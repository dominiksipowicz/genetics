#!/usr/bin/env python3
"""
Compare CRAM Files

This script compares CRAM files to analyze:
- Shared genomic regions with coverage
- Variant calling and comparison
- DNA sharing analysis between samples
- Inheritance patterns

Requires samtools, bcftools, and optionally PLINK for kinship analysis.
"""

import subprocess
import pandas as pd
import numpy as np
import json
import sys
import os
from pathlib import Path
from collections import defaultdict
import tempfile


class CramComparator:
    def __init__(self, cram1_path, cram2_path, sample1_name, sample2_name, reference_path=None):
        """
        Initialize CRAM comparator
        
        Args:
            cram1_path (str): Path to first CRAM file
            cram2_path (str): Path to second CRAM file
            sample1_name (str): Name for first sample
            sample2_name (str): Name for second sample
            reference_path (str): Path to reference genome
        """
        self.cram1_path = Path(cram1_path)
        self.cram2_path = Path(cram2_path)
        self.sample1_name = sample1_name
        self.sample2_name = sample2_name
        self.reference_path = Path(reference_path) if reference_path else None
        
        # Validate files exist
        for path, name in [(self.cram1_path, sample1_name), (self.cram2_path, sample2_name)]:
            if not path.exists():
                raise FileNotFoundError(f"CRAM file not found for {name}: {path}")
    
    def check_tools(self):
        """Check if required tools are available"""
        tools = {
            'samtools': ['samtools', '--version'],
            'bcftools': ['bcftools', '--version']
        }
        
        available = {}
        for tool, cmd in tools.items():
            try:
                result = subprocess.run(cmd, capture_output=True, text=True)
                available[tool] = result.returncode == 0
            except FileNotFoundError:
                available[tool] = False
        
        return available
    
    def get_shared_coverage_regions(self, min_depth=10):
        """
        Find genomic regions where both samples have coverage above threshold
        """
        print(f"Finding shared coverage regions (min depth: {min_depth}x)...")
        
        shared_regions = []
        
        try:
            # Get coverage for both samples on main chromosomes
            chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
            
            for chrom in chromosomes:
                print(f"  Analyzing chromosome {chrom}...")
                
                # Get coverage for sample 1
                cmd1 = ['samtools', 'depth', '-r', f'{chrom}', str(self.cram1_path)]
                if self.reference_path:
                    cmd1.extend(['-f', str(self.reference_path)])
                
                # Get coverage for sample 2
                cmd2 = ['samtools', 'depth', '-r', f'{chrom}', str(self.cram2_path)]
                if self.reference_path:
                    cmd2.extend(['-f', str(self.reference_path)])
                
                # For demonstration, we'll sample positions to avoid memory issues
                # In production, you might want to process in chunks
                try:
                    # Sample first 10000 positions per chromosome
                    result1 = subprocess.run(cmd1 + ['|', 'head', '-10000'], 
                                           shell=True, capture_output=True, text=True, timeout=60)
                    result2 = subprocess.run(cmd2 + ['|', 'head', '-10000'], 
                                           shell=True, capture_output=True, text=True, timeout=60)
                    
                    if result1.returncode == 0 and result2.returncode == 0:
                        # Parse coverage data
                        coverage1 = {}
                        for line in result1.stdout.strip().split('\n'):
                            if line:
                                parts = line.split('\t')
                                if len(parts) >= 3:
                                    pos = int(parts[1])
                                    depth = int(parts[2])
                                    coverage1[pos] = depth
                        
                        coverage2 = {}
                        for line in result2.stdout.strip().split('\n'):
                            if line:
                                parts = line.split('\t')
                                if len(parts) >= 3:
                                    pos = int(parts[1])
                                    depth = int(parts[2])
                                    coverage2[pos] = depth
                        
                        # Find shared positions with sufficient coverage
                        shared_positions = []
                        for pos in coverage1:
                            if pos in coverage2:
                                if coverage1[pos] >= min_depth and coverage2[pos] >= min_depth:
                                    shared_positions.append({
                                        'chromosome': chrom,
                                        'position': pos,
                                        f'{self.sample1_name}_depth': coverage1[pos],
                                        f'{self.sample2_name}_depth': coverage2[pos]
                                    })
                        
                        shared_regions.extend(shared_positions)
                        print(f"    Found {len(shared_positions)} shared positions")
                
                except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
                    print(f"    Error processing chromosome {chrom}: {e}")
                    continue
        
        except Exception as e:
            print(f"Error in coverage analysis: {e}")
        
        return pd.DataFrame(shared_regions)
    
    def call_variants(self, output_dir="variant_calling"):
        """
        Call variants for both samples and compare
        """
        print("Calling variants for comparison...")
        
        os.makedirs(output_dir, exist_ok=True)
        
        if not self.reference_path or not self.reference_path.exists():
            print("Warning: Reference genome required for variant calling")
            return None
        
        variant_files = {}
        
        for sample_name, cram_path in [(self.sample1_name, self.cram1_path), 
                                       (self.sample2_name, self.cram2_path)]:
            print(f"  Calling variants for {sample_name}...")
            
            vcf_path = Path(output_dir) / f"{sample_name}_variants.vcf"
            
            try:
                # Call variants using bcftools mpileup + call
                cmd = [
                    'bcftools', 'mpileup',
                    '-f', str(self.reference_path),
                    '-Ou',  # Uncompressed BCF output
                    str(cram_path)
                ]
                
                # Pipe to bcftools call
                cmd_call = [
                    'bcftools', 'call',
                    '-mv',  # Multiallelic and variants only
                    '-Ov',  # VCF output
                    '-o', str(vcf_path)
                ]
                
                # Run mpileup | call pipeline
                process1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                process2 = subprocess.Popen(cmd_call, stdin=process1.stdout, 
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process1.stdout.close()
                
                stdout, stderr = process2.communicate(timeout=300)  # 5 minute timeout
                
                if process2.returncode == 0 and vcf_path.exists():
                    variant_files[sample_name] = vcf_path
                    print(f"    Variants saved to: {vcf_path}")
                else:
                    print(f"    Error calling variants for {sample_name}: {stderr.decode()}")
            
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError) as e:
                print(f"    Error calling variants for {sample_name}: {e}")
        
        return variant_files
    
    def compare_variants(self, variant_files):
        """
        Compare variants between samples
        """
        if len(variant_files) != 2:
            print("Need variants for both samples to compare")
            return None
        
        print("Comparing variants between samples...")
        
        variants = {}
        
        # Parse VCF files
        for sample_name, vcf_path in variant_files.items():
            print(f"  Parsing variants for {sample_name}...")
            sample_variants = []
            
            try:
                with open(vcf_path, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        
                        parts = line.strip().split('\t')
                        if len(parts) >= 8:
                            chrom = parts[0]
                            pos = int(parts[1])
                            ref = parts[3]
                            alt = parts[4]
                            qual = parts[5]
                            
                            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
                            sample_variants.append({
                                'id': variant_id,
                                'chromosome': chrom,
                                'position': pos,
                                'reference': ref,
                                'alternate': alt,
                                'quality': qual
                            })
                
                variants[sample_name] = sample_variants
                print(f"    Found {len(sample_variants)} variants")
            
            except Exception as e:
                print(f"    Error parsing variants for {sample_name}: {e}")
        
        # Compare variants
        if len(variants) == 2:
            sample_names = list(variants.keys())
            variants1 = {v['id'] for v in variants[sample_names[0]]}
            variants2 = {v['id'] for v in variants[sample_names[1]]}
            
            shared_variants = variants1 & variants2
            unique_to_1 = variants1 - variants2
            unique_to_2 = variants2 - variants1
            
            comparison = {
                'total_variants': {
                    sample_names[0]: len(variants1),
                    sample_names[1]: len(variants2)
                },
                'shared_variants': len(shared_variants),
                'unique_variants': {
                    sample_names[0]: len(unique_to_1),
                    sample_names[1]: len(unique_to_2)
                },
                'jaccard_similarity': len(shared_variants) / len(variants1 | variants2) if variants1 | variants2 else 0
            }
            
            return comparison
        
        return None
    
    def estimate_dna_sharing(self, shared_coverage_df, variant_comparison):
        """
        Estimate DNA sharing between samples
        """
        print("Estimating DNA sharing...")
        
        sharing_metrics = {}
        
        # Coverage-based sharing
        if not shared_coverage_df.empty:
            total_positions = len(shared_coverage_df)
            
            # Calculate correlation in coverage depth
            depths1 = shared_coverage_df[f'{self.sample1_name}_depth']
            depths2 = shared_coverage_df[f'{self.sample2_name}_depth']
            
            if len(depths1) > 1:
                correlation = np.corrcoef(depths1, depths2)[0, 1]
                sharing_metrics['coverage_correlation'] = correlation
            
            sharing_metrics['shared_positions_analyzed'] = total_positions
        
        # Variant-based sharing
        if variant_comparison:
            sharing_metrics.update(variant_comparison)
        
        # Estimate relationship
        if 'jaccard_similarity' in sharing_metrics:
            jaccard = sharing_metrics['jaccard_similarity']
            if jaccard > 0.4:
                relationship = "Close relative (parent-child or sibling)"
            elif jaccard > 0.2:
                relationship = "Distant relative"
            else:
                relationship = "Unrelated or very distant"
            
            sharing_metrics['estimated_relationship'] = relationship
        
        return sharing_metrics
    
    def run_comparison(self):
        """
        Run complete comparison analysis
        """
        print(f"\n=== Comparing CRAM files ===")
        print(f"{self.sample1_name}: {self.cram1_path}")
        print(f"{self.sample2_name}: {self.cram2_path}")
        
        # Check tools
        tools = self.check_tools()
        missing_tools = [tool for tool, available in tools.items() if not available]
        
        if missing_tools:
            print(f"Warning: Missing tools: {missing_tools}")
            print("Some analyses may not be available.")
        
        results = {
            'sample1': self.sample1_name,
            'sample2': self.sample2_name,
            'cram1_path': str(self.cram1_path),
            'cram2_path': str(self.cram2_path)
        }
        
        # 1. Find shared coverage regions
        try:
            shared_coverage = self.get_shared_coverage_regions()
            results['shared_coverage'] = {
                'total_positions': len(shared_coverage),
                'chromosomes_analyzed': shared_coverage['chromosome'].nunique() if not shared_coverage.empty else 0
            }
            
            if not shared_coverage.empty:
                # Save shared positions
                coverage_file = f"{self.sample1_name}_{self.sample2_name}_shared_coverage.csv"
                shared_coverage.to_csv(coverage_file, index=False)
                print(f"Shared coverage positions saved to: {coverage_file}")
        
        except Exception as e:
            print(f"Error in coverage analysis: {e}")
            shared_coverage = pd.DataFrame()
            results['shared_coverage'] = {'error': str(e)}
        
        # 2. Call and compare variants
        variant_comparison = None
        if tools.get('bcftools', False):
            try:
                variant_files = self.call_variants()
                if variant_files:
                    variant_comparison = self.compare_variants(variant_files)
                    results['variant_comparison'] = variant_comparison
            except Exception as e:
                print(f"Error in variant analysis: {e}")
                results['variant_comparison'] = {'error': str(e)}
        else:
            print("Skipping variant analysis (bcftools not available)")
        
        # 3. Estimate DNA sharing
        try:
            sharing_metrics = self.estimate_dna_sharing(shared_coverage, variant_comparison)
            results['dna_sharing'] = sharing_metrics
        except Exception as e:
            print(f"Error in DNA sharing analysis: {e}")
            results['dna_sharing'] = {'error': str(e)}
        
        return results
    
    def print_summary(self, results):
        """Print human-readable summary"""
        print(f"\n=== Comparison Summary ===")
        print(f"Samples: {results['sample1']} vs {results['sample2']}")
        
        # Coverage summary
        if 'shared_coverage' in results and 'total_positions' in results['shared_coverage']:
            coverage = results['shared_coverage']
            print(f"\nShared Coverage:")
            print(f"  Positions with coverage in both samples: {coverage['total_positions']:,}")
            print(f"  Chromosomes analyzed: {coverage['chromosomes_analyzed']}")
        
        # Variant summary
        if 'variant_comparison' in results and 'shared_variants' in results['variant_comparison']:
            variants = results['variant_comparison']
            print(f"\nVariant Comparison:")
            print(f"  Shared variants: {variants['shared_variants']:,}")
            print(f"  Jaccard similarity: {variants['jaccard_similarity']:.3f}")
            
            for sample, count in variants['total_variants'].items():
                print(f"  Total variants in {sample}: {count:,}")
        
        # DNA sharing summary
        if 'dna_sharing' in results:
            sharing = results['dna_sharing']
            print(f"\nDNA Sharing Analysis:")
            
            if 'coverage_correlation' in sharing:
                print(f"  Coverage correlation: {sharing['coverage_correlation']:.3f}")
            
            if 'estimated_relationship' in sharing:
                print(f"  Estimated relationship: {sharing['estimated_relationship']}")


def main():
    """Main function"""
    
    # Define sample paths
    dominik_cram = '/Users/dom/work/genetics/Dom/Dominik-NG173LPFBH.cram'
    jadwiga_cram = '/Users/dom/work/genetics/Jadwiga/MG46K47NL.cram'
    reference_path = '/Users/dom/work/genetics/Dom/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    
    # Check if files exist
    if not os.path.exists(dominik_cram):
        print(f"Error: Dominik's CRAM file not found: {dominik_cram}")
        return
    
    if not os.path.exists(jadwiga_cram):
        print(f"Error: Jadwiga's CRAM file not found: {jadwiga_cram}")
        return
    
    if not os.path.exists(reference_path):
        print(f"Warning: Reference genome not found: {reference_path}")
        reference_path = None
    
    # Run comparison
    comparator = CramComparator(
        dominik_cram, jadwiga_cram,
        'Dominik', 'Jadwiga',
        reference_path
    )
    
    results = comparator.run_comparison()
    comparator.print_summary(results)
    
    # Save results
    output_file = 'dominik_jadwiga_comparison.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\nDetailed results saved to: {output_file}")
    print("\n=== Analysis Complete ===")


if __name__ == "__main__":
    main()



