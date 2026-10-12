#!/usr/bin/env python3
"""
Direct CRAM Analysis - Single Source of Truth

This script extracts variants directly from CRAM files without relying on 
pre-existing VCF data. It performs:
1. Variant calling directly from CRAM files
2. Gene annotation of variants
3. Comparison of genes between samples

Requires: samtools, bcftools, and optionally VEP or SnpEff for annotation
"""

import subprocess
import pandas as pd
import json
import os
import sys
from pathlib import Path
import tempfile
import shutil


class DirectCramAnalyzer:
    def __init__(self, reference_path):
        """
        Initialize with reference genome path
        
        Args:
            reference_path (str): Path to reference genome FASTA file
        """
        self.reference_path = Path(reference_path)
        if not self.reference_path.exists():
            raise FileNotFoundError(f"Reference genome not found: {reference_path}")
        
        self.samples = {}
        self.output_dir = Path("cram_direct_output")
        self.output_dir.mkdir(exist_ok=True)
        
    def add_sample(self, name, cram_path):
        """Add a sample for analysis"""
        cram_path = Path(cram_path)
        if not cram_path.exists():
            raise FileNotFoundError(f"CRAM file not found: {cram_path}")
        
        self.samples[name] = cram_path
        print(f"Added sample: {name} -> {cram_path}")
    
    def check_tools(self):
        """Check if required tools are available"""
        tools = ['samtools', 'bcftools']
        available = {}
        
        for tool in tools:
            try:
                result = subprocess.run([tool, '--version'], 
                                      capture_output=True, text=True)
                available[tool] = result.returncode == 0
                if available[tool]:
                    version = result.stdout.split('\n')[0]
                    print(f"✓ {tool}: {version}")
                else:
                    print(f"✗ {tool}: Not working properly")
            except FileNotFoundError:
                available[tool] = False
                print(f"✗ {tool}: Not found")
        
        return available
    
    def call_variants_from_cram(self, sample_name, chromosomes=None, max_depth=250):
        """
        Call variants directly from CRAM file
        
        Args:
            sample_name (str): Name of the sample
            chromosomes (list): List of chromosomes to analyze (default: 1-22, X, Y)
            max_depth (int): Maximum read depth for mpileup
        """
        if sample_name not in self.samples:
            raise ValueError(f"Sample {sample_name} not found")
        
        cram_path = self.samples[sample_name]
        output_vcf = self.output_dir / f"{sample_name}_variants.vcf.gz"
        
        print(f"\nCalling variants for {sample_name}...")
        print(f"  CRAM: {cram_path}")
        print(f"  Output: {output_vcf}")
        
        # Default chromosomes if not specified
        if chromosomes is None:
            chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
        
        try:
            # Create chromosome regions string
            regions = ' '.join([f"chr{chrom}" for chrom in chromosomes])
            
            # Build the command
            # Use mpileup with specific regions to avoid memory issues
            cmd = [
                'bcftools', 'mpileup',
                '-f', str(self.reference_path),
                '-r', regions,  # Specify regions
                '-d', str(max_depth),  # Max depth
                '-q', '20',  # Min mapping quality
                '-Q', '20',  # Min base quality
                '-a', 'FORMAT/AD,FORMAT/DP',  # Additional format fields
                '-Ou',  # Uncompressed BCF output
                str(cram_path)
            ]
            
            # Pipe to bcftools call
            call_cmd = [
                'bcftools', 'call',
                '-mv',  # Multiallelic and variants only
                '-Oz',  # Compressed VCF output
                '-o', str(output_vcf)
            ]
            
            print(f"  Running: {' '.join(cmd)} | {' '.join(call_cmd)}")
            
            # Execute the pipeline
            with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc1:
                with subprocess.Popen(call_cmd, stdin=proc1.stdout, 
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc2:
                    proc1.stdout.close()
                    stdout, stderr = proc2.communicate(timeout=1800)  # 30 minute timeout
                    
                    if proc2.returncode != 0:
                        print(f"  Error in variant calling: {stderr.decode()}")
                        return None
            
            # Index the VCF file
            subprocess.run(['bcftools', 'index', str(output_vcf)], check=True)
            
            print(f"  ✓ Variants called successfully")
            return output_vcf
            
        except subprocess.TimeoutExpired:
            print(f"  ✗ Variant calling timed out for {sample_name}")
            return None
        except subprocess.CalledProcessError as e:
            print(f"  ✗ Error calling variants for {sample_name}: {e}")
            return None
    
    def parse_vcf_file(self, vcf_path):
        """
        Parse VCF file and extract variant information
        
        Args:
            vcf_path (Path): Path to VCF file
            
        Returns:
            pd.DataFrame: DataFrame with variant information
        """
        print(f"Parsing VCF file: {vcf_path}")
        
        variants = []
        
        try:
            # Use bcftools query to extract information
            cmd = [
                'bcftools', 'query',
                '-f', '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\n',
                str(vcf_path)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        variants.append({
                            'CHROM': parts[0],
                            'POS': int(parts[1]),
                            'REF': parts[2],
                            'ALT': parts[3],
                            'QUAL': float(parts[4]) if parts[4] != '.' else None,
                            'DP': int(parts[5]) if len(parts) > 5 and parts[5] != '.' else None
                        })
            
            df = pd.DataFrame(variants)
            
            if not df.empty:
                # Add variant type
                df['variant_type'] = df.apply(
                    lambda row: 'SNP' if len(row['REF']) == 1 and len(row['ALT']) == 1 else 'INDEL', 
                    axis=1
                )
                
                # Create variant ID
                df['variant_id'] = df['CHROM'] + ':' + df['POS'].astype(str) + ':' + df['REF'] + ':' + df['ALT']
            
            print(f"  Parsed {len(df)} variants")
            return df
            
        except subprocess.CalledProcessError as e:
            print(f"  Error parsing VCF: {e}")
            return pd.DataFrame()
    
    def annotate_variants_with_genes(self, variants_df, sample_name):
        """
        Annotate variants with gene information using a simple approach
        This is a basic implementation - in production, you'd use VEP or SnpEff
        """
        print(f"Annotating variants for {sample_name}...")
        
        # For this demonstration, we'll create a simple gene annotation
        # In practice, you would use tools like VEP, SnpEff, or ANNOVAR
        
        if variants_df.empty:
            return variants_df
        
        # Simple gene region approximation (this is very basic)
        # In reality, you'd use a proper gene annotation database
        annotated_variants = variants_df.copy()
        
        # Mock gene annotation based on chromosome and position
        # This is just for demonstration - replace with real annotation
        def simple_gene_annotation(chrom, pos):
            """Very basic gene region estimation"""
            # This is a simplified mock - in reality use proper gene databases
            gene_regions = {
                'chr1': [(1000000, 2000000, 'GENE1A'), (5000000, 6000000, 'GENE1B')],
                'chr2': [(1000000, 1500000, 'GENE2A'), (3000000, 4000000, 'GENE2B')],
                # Add more as needed
            }
            
            if chrom in gene_regions:
                for start, end, gene in gene_regions[chrom]:
                    if start <= pos <= end:
                        return gene
            
            return f"INTERGENIC_{chrom}_{pos//1000000}MB"
        
        annotated_variants['gene'] = annotated_variants.apply(
            lambda row: simple_gene_annotation(row['CHROM'], row['POS']), axis=1
        )
        
        # Count genes
        gene_counts = annotated_variants['gene'].value_counts()
        print(f"  Found variants in {len(gene_counts)} gene regions")
        
        return annotated_variants
    
    def compare_genes(self, sample1_variants, sample2_variants, sample1_name, sample2_name):
        """Compare genes between two samples"""
        print(f"\nComparing genes between {sample1_name} and {sample2_name}...")
        
        if sample1_variants.empty or sample2_variants.empty:
            print("Cannot compare - one or both samples have no variants")
            return {}
        
        # Get genes for each sample
        sample1_genes = set(sample1_variants['gene'].unique())
        sample2_genes = set(sample2_variants['gene'].unique())
        
        # Find overlaps
        shared_genes = sample1_genes & sample2_genes
        sample1_unique = sample1_genes - sample2_genes
        sample2_unique = sample2_genes - sample1_genes
        
        comparison = {
            f'{sample1_name}_genes': len(sample1_genes),
            f'{sample2_name}_genes': len(sample2_genes),
            'shared_genes': len(shared_genes),
            f'{sample1_name}_unique_genes': len(sample1_unique),
            f'{sample2_name}_unique_genes': len(sample2_unique),
            'shared_gene_list': sorted(list(shared_genes)),
            f'{sample1_name}_unique_gene_list': sorted(list(sample1_unique)),
            f'{sample2_name}_unique_gene_list': sorted(list(sample2_unique))
        }
        
        print(f"  {sample1_name} genes: {comparison[f'{sample1_name}_genes']}")
        print(f"  {sample2_name} genes: {comparison[f'{sample2_name}_genes']}")
        print(f"  Shared genes: {comparison['shared_genes']}")
        print(f"  {sample1_name} unique: {comparison[f'{sample1_name}_unique_genes']}")
        print(f"  {sample2_name} unique: {comparison[f'{sample2_name}_unique_genes']}")
        
        return comparison
    
    def analyze_all_samples(self, chromosomes=None):
        """Run complete analysis on all samples"""
        print("="*60)
        print("DIRECT CRAM ANALYSIS - SINGLE SOURCE OF TRUTH")
        print("="*60)
        
        # Check tools
        tools = self.check_tools()
        if not all(tools.values()):
            print("Missing required tools. Please install samtools and bcftools.")
            return None
        
        results = {
            'reference_genome': str(self.reference_path),
            'samples': {},
            'comparisons': {}
        }
        
        # Process each sample
        sample_variants = {}
        
        for sample_name in self.samples:
            print(f"\n{'='*40}")
            print(f"Processing {sample_name}")
            print(f"{'='*40}")
            
            # Call variants
            vcf_file = self.call_variants_from_cram(sample_name, chromosomes)
            
            if vcf_file and vcf_file.exists():
                # Parse variants
                variants_df = self.parse_vcf_file(vcf_file)
                
                if not variants_df.empty:
                    # Annotate with genes
                    annotated_variants = self.annotate_variants_with_genes(variants_df, sample_name)
                    sample_variants[sample_name] = annotated_variants
                    
                    # Store results
                    results['samples'][sample_name] = {
                        'vcf_file': str(vcf_file),
                        'total_variants': len(annotated_variants),
                        'snps': len(annotated_variants[annotated_variants['variant_type'] == 'SNP']),
                        'indels': len(annotated_variants[annotated_variants['variant_type'] == 'INDEL']),
                        'genes_affected': len(annotated_variants['gene'].unique()),
                        'chromosomes': len(annotated_variants['CHROM'].unique())
                    }
                    
                    # Save detailed variants
                    variant_file = self.output_dir / f"{sample_name}_annotated_variants.csv"
                    annotated_variants.to_csv(variant_file, index=False)
                    print(f"  Saved annotated variants to: {variant_file}")
                else:
                    print(f"  No variants found for {sample_name}")
            else:
                print(f"  Failed to call variants for {sample_name}")
        
        # Compare between samples
        sample_names = list(sample_variants.keys())
        if len(sample_names) >= 2:
            for i in range(len(sample_names)):
                for j in range(i + 1, len(sample_names)):
                    sample1, sample2 = sample_names[i], sample_names[j]
                    comparison = self.compare_genes(
                        sample_variants[sample1], 
                        sample_variants[sample2],
                        sample1, sample2
                    )
                    results['comparisons'][f"{sample1}_vs_{sample2}"] = comparison
        
        # Save results
        results_file = self.output_dir / "cram_direct_analysis_results.json"
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        print(f"\n{'='*60}")
        print("ANALYSIS COMPLETE")
        print(f"{'='*60}")
        print(f"Results saved to: {results_file}")
        
        return results


def main():
    """Main function"""
    
    # Configuration
    reference_genome = "/Users/dom/work/genetics/Dom/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    
    samples = {
        'Dominik': '/Users/dom/work/genetics/Dom/Dominik-NG173LPFBH.cram',
        'Jadwiga': '/Users/dom/work/genetics/Jadwiga/MG46K47NL.cram'
    }
    
    # Limit to main chromosomes for faster analysis
    chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                   '11', '12', '13', '14', '15', '16', '17', '18', '19', 
                   '20', '21', '22', 'X', 'Y']
    
    try:
        # Initialize analyzer
        analyzer = DirectCramAnalyzer(reference_genome)
        
        # Add samples
        for name, path in samples.items():
            if os.path.exists(path):
                analyzer.add_sample(name, path)
            else:
                print(f"Warning: CRAM file not found for {name}: {path}")
        
        # Run analysis
        if analyzer.samples:
            results = analyzer.analyze_all_samples(chromosomes)
            
            # Print summary
            if results and 'samples' in results:
                print(f"\n{'='*60}")
                print("SUMMARY")
                print(f"{'='*60}")
                
                for sample_name, stats in results['samples'].items():
                    print(f"\n{sample_name}:")
                    print(f"  Total variants: {stats['total_variants']:,}")
                    print(f"  SNPs: {stats['snps']:,}")
                    print(f"  INDELs: {stats['indels']:,}")
                    print(f"  Genes affected: {stats['genes_affected']:,}")
                    print(f"  Chromosomes: {stats['chromosomes']}")
                
                if 'comparisons' in results:
                    print(f"\nComparisons:")
                    for comp_name, comp_data in results['comparisons'].items():
                        print(f"\n{comp_name}:")
                        print(f"  Shared genes: {comp_data['shared_genes']}")
                        print(f"  Total unique genes: {comp_data.get('Dominik_genes', 0) + comp_data.get('Jadwiga_genes', 0) - comp_data['shared_genes']}")
        else:
            print("No valid CRAM files found for analysis")
            
    except Exception as e:
        print(f"Error in analysis: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()



