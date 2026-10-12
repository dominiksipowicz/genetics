#!/usr/bin/env python3
"""
CRAM Statistics Analysis

This script analyzes CRAM files to extract statistics about:
- Number of reads and positions
- Chromosome coverage
- Mapping quality statistics
- Coverage depth statistics

Requires samtools to be installed and available in PATH.
"""

import subprocess
import pandas as pd
import json
import sys
import os
from pathlib import Path


class CramAnalyzer:
    def __init__(self, cram_path, reference_path=None):
        """
        Initialize CRAM analyzer
        
        Args:
            cram_path (str): Path to CRAM file
            reference_path (str): Path to reference genome (optional)
        """
        self.cram_path = Path(cram_path)
        self.reference_path = Path(reference_path) if reference_path else None
        self.sample_name = self.cram_path.stem
        
        if not self.cram_path.exists():
            raise FileNotFoundError(f"CRAM file not found: {cram_path}")
    
    def check_samtools(self):
        """Check if samtools is available"""
        try:
            result = subprocess.run(['samtools', '--version'], 
                                  capture_output=True, text=True)
            return result.returncode == 0
        except FileNotFoundError:
            return False
    
    def get_basic_stats(self):
        """Get basic statistics from CRAM file using samtools stats"""
        print(f"Analyzing basic statistics for {self.sample_name}...")
        
        cmd = ['samtools', 'stats', str(self.cram_path)]
        if self.reference_path:
            cmd.extend(['--reference', str(self.reference_path)])
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            stats = {}
            
            for line in result.stdout.split('\n'):
                if line.startswith('SN'):
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        key = parts[1].strip(':')
                        try:
                            value = int(parts[2])
                        except ValueError:
                            try:
                                value = float(parts[2])
                            except ValueError:
                                value = parts[2]
                        stats[key] = value
            
            return stats
        except subprocess.CalledProcessError as e:
            print(f"Error running samtools stats: {e}")
            return {}
    
    def get_chromosome_coverage(self):
        """Get coverage statistics per chromosome"""
        print(f"Analyzing chromosome coverage for {self.sample_name}...")
        
        cmd = ['samtools', 'idxstats', str(self.cram_path)]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            chromosomes = []
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        chrom = parts[0]
                        length = int(parts[1])
                        mapped_reads = int(parts[2])
                        unmapped_reads = int(parts[3])
                        
                        chromosomes.append({
                            'chromosome': chrom,
                            'length': length,
                            'mapped_reads': mapped_reads,
                            'unmapped_reads': unmapped_reads,
                            'total_reads': mapped_reads + unmapped_reads
                        })
            
            return pd.DataFrame(chromosomes)
        except subprocess.CalledProcessError as e:
            print(f"Error running samtools idxstats: {e}")
            return pd.DataFrame()
    
    def get_coverage_depth(self, regions=None):
        """Get coverage depth statistics"""
        print(f"Analyzing coverage depth for {self.sample_name}...")
        
        cmd = ['samtools', 'depth', str(self.cram_path)]
        if regions:
            cmd.extend(['-r', regions])
        
        try:
            # For large files, we'll sample to avoid memory issues
            result = subprocess.run(cmd + ['|', 'head', '-100000'], 
                                  shell=True, capture_output=True, text=True, check=True)
            
            depths = []
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        depths.append(int(parts[2]))
            
            if depths:
                return {
                    'mean_depth': sum(depths) / len(depths),
                    'max_depth': max(depths),
                    'min_depth': min(depths),
                    'positions_sampled': len(depths)
                }
            return {}
        except subprocess.CalledProcessError as e:
            print(f"Error running samtools depth: {e}")
            return {}
    
    def analyze(self):
        """Run complete analysis"""
        if not self.check_samtools():
            print("Error: samtools not found. Please install samtools to analyze CRAM files.")
            return None
        
        print(f"\n=== CRAM Analysis for {self.sample_name} ===")
        print(f"File: {self.cram_path}")
        
        # Basic statistics
        basic_stats = self.get_basic_stats()
        
        # Chromosome coverage
        chrom_df = self.get_chromosome_coverage()
        
        # Coverage depth (sample)
        depth_stats = self.get_coverage_depth()
        
        # Compile results
        results = {
            'sample_name': self.sample_name,
            'cram_path': str(self.cram_path),
            'basic_stats': basic_stats,
            'chromosome_stats': chrom_df.to_dict('records') if not chrom_df.empty else [],
            'depth_stats': depth_stats
        }
        
        return results
    
    def print_summary(self, results):
        """Print a human-readable summary"""
        if not results:
            return
        
        print(f"\n=== Summary for {results['sample_name']} ===")
        
        # Basic stats
        basic = results['basic_stats']
        if basic:
            print(f"Total reads: {basic.get('raw total sequences', 'N/A'):,}")
            print(f"Mapped reads: {basic.get('reads mapped', 'N/A'):,}")
            print(f"Mapping rate: {basic.get('percentage of properly paired reads (%)', 'N/A')}%")
            print(f"Average read length: {basic.get('average length', 'N/A')}")
            print(f"Total bases: {basic.get('total length', 'N/A'):,}")
        
        # Chromosome stats
        chrom_stats = results['chromosome_stats']
        if chrom_stats:
            df = pd.DataFrame(chrom_stats)
            # Filter out unmapped and small contigs
            main_chroms = df[df['chromosome'].str.match(r'^(chr)?([0-9]+|X|Y|MT?)$', na=False)]
            
            print(f"\nChromosome Coverage:")
            print(f"Total chromosomes/contigs: {len(df)}")
            print(f"Main chromosomes (1-22, X, Y, MT): {len(main_chroms)}")
            
            if not main_chroms.empty:
                print(f"Total mapped reads on main chromosomes: {main_chroms['mapped_reads'].sum():,}")
                print(f"Most covered chromosome: {main_chroms.loc[main_chroms['mapped_reads'].idxmax(), 'chromosome']}")
        
        # Depth stats
        depth = results['depth_stats']
        if depth:
            print(f"\nCoverage Depth (sampled):")
            print(f"Mean depth: {depth['mean_depth']:.2f}x")
            print(f"Max depth: {depth['max_depth']}x")
            print(f"Positions sampled: {depth['positions_sampled']:,}")


def main():
    """Main function to run CRAM analysis"""
    
    # Define CRAM file paths
    cram_files = {
        'Dominik': '/Users/dom/work/genetics/Dom/Dominik-NG173LPFBH.cram',
        'Jadwiga': '/Users/dom/work/genetics/Jadwiga/MG46K47NL.cram'
    }
    
    # Reference genome path (if available)
    reference_path = '/Users/dom/work/genetics/Dom/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    if not os.path.exists(reference_path):
        reference_path = None
        print("Warning: Reference genome not found. Analysis may be limited.")
    
    results = {}
    
    # Analyze each CRAM file
    for sample_name, cram_path in cram_files.items():
        if os.path.exists(cram_path):
            print(f"\n{'='*50}")
            print(f"Analyzing {sample_name}")
            print(f"{'='*50}")
            
            analyzer = CramAnalyzer(cram_path, reference_path)
            result = analyzer.analyze()
            
            if result:
                analyzer.print_summary(result)
                results[sample_name] = result
                
                # Save detailed results to JSON
                output_file = f"{sample_name}_cram_stats.json"
                with open(output_file, 'w') as f:
                    json.dump(result, f, indent=2, default=str)
                print(f"\nDetailed results saved to: {output_file}")
        else:
            print(f"CRAM file not found for {sample_name}: {cram_path}")
    
    # Save combined results
    if results:
        with open('all_cram_stats.json', 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\nCombined results saved to: all_cram_stats.json")
    
    print(f"\n{'='*50}")
    print("Analysis complete!")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()
