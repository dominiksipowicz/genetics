#!/usr/bin/env python3
"""
Simple CRAM Analysis and Comparison

A more practical approach to analyzing CRAM files that focuses on:
1. Basic statistics from existing data
2. Using pre-computed VCF files for variant comparison
3. Estimating DNA sharing using available genetic data

This script is designed to work with the existing VCF files and provide
meaningful insights about genetic relatedness.
"""

import pandas as pd
import json
import os
from pathlib import Path
import subprocess


def analyze_existing_vcf_data():
    """Analyze the existing parsed VCF data for both samples"""
    
    # Paths to existing VCF data
    dominik_vcf_csv = "/Users/dom/work/genetics/Dom/parsed_vcf_data.csv"
    jadwiga_vcf_csv = "/Users/dom/work/genetics/Jadwiga/parsed_vcf_data.csv"
    
    results = {}
    
    # Load and analyze each sample's VCF data
    for sample_name, vcf_path in [("Dominik", dominik_vcf_csv), ("Jadwiga", jadwiga_vcf_csv)]:
        print(f"\nAnalyzing {sample_name}'s VCF data...")
        
        if not os.path.exists(vcf_path):
            print(f"  Warning: VCF data not found for {sample_name}")
            continue
        
        try:
            df = pd.read_csv(vcf_path)
            
            # Basic statistics
            stats = {
                'total_variants': len(df),
                'chromosomes': df['CHROM'].nunique() if 'CHROM' in df.columns else 'N/A',
                'unique_positions': df.groupby(['CHROM', 'POS']).size().shape[0] if all(col in df.columns for col in ['CHROM', 'POS']) else 'N/A',
                'variant_types': {},
                'quality_stats': {}
            }
            
            # Analyze variant types if REF and ALT columns exist
            if all(col in df.columns for col in ['REF', 'ALT']):
                # SNPs vs INDELs
                df['variant_type'] = df.apply(lambda row: 'SNP' if len(row['REF']) == 1 and len(row['ALT']) == 1 else 'INDEL', axis=1)
                stats['variant_types'] = df['variant_type'].value_counts().to_dict()
            
            # Quality statistics if available
            if 'QUAL' in df.columns:
                stats['quality_stats'] = {
                    'mean_quality': df['QUAL'].mean(),
                    'median_quality': df['QUAL'].median(),
                    'min_quality': df['QUAL'].min(),
                    'max_quality': df['QUAL'].max()
                }
            
            # Chromosome distribution
            if 'CHROM' in df.columns:
                chrom_counts = df['CHROM'].value_counts()
                main_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
                main_chrom_counts = {chrom: chrom_counts.get(chrom, 0) for chrom in main_chroms if chrom in chrom_counts}
                stats['chromosome_distribution'] = dict(sorted(main_chrom_counts.items(), key=lambda x: x[1], reverse=True)[:10])
            
            results[sample_name] = {
                'vcf_file': vcf_path,
                'statistics': stats,
                'sample_data': df.head().to_dict() if len(df) > 0 else {}
            }
            
            print(f"  Total variants: {stats['total_variants']:,}")
            print(f"  Chromosomes: {stats['chromosomes']}")
            if stats['variant_types']:
                for vtype, count in stats['variant_types'].items():
                    print(f"  {vtype}s: {count:,}")
        
        except Exception as e:
            print(f"  Error analyzing {sample_name}: {e}")
            results[sample_name] = {'error': str(e)}
    
    return results


def compare_variants(dominik_data, jadwiga_data):
    """Compare variants between Dominik and Jadwiga using VCF data"""
    
    print("\n=== Comparing Variants ===")
    
    if 'error' in dominik_data or 'error' in jadwiga_data:
        print("Cannot compare variants due to data loading errors")
        return {}
    
    try:
        # Load the actual dataframes for comparison
        dominik_df = pd.read_csv(dominik_data['vcf_file'])
        jadwiga_df = pd.read_csv(jadwiga_data['vcf_file'])
        
        # Create variant identifiers
        def create_variant_id(df):
            if all(col in df.columns for col in ['CHROM', 'POS', 'REF', 'ALT']):
                return df['CHROM'].astype(str) + ':' + df['POS'].astype(str) + ':' + df['REF'] + ':' + df['ALT']
            else:
                return pd.Series([])
        
        dominik_variants = set(create_variant_id(dominik_df))
        jadwiga_variants = set(create_variant_id(jadwiga_df))
        
        # Remove empty identifiers
        dominik_variants.discard('')
        jadwiga_variants.discard('')
        
        if not dominik_variants or not jadwiga_variants:
            print("Unable to create variant identifiers - missing required columns")
            return {}
        
        # Calculate overlap
        shared_variants = dominik_variants & jadwiga_variants
        dominik_unique = dominik_variants - jadwiga_variants
        jadwiga_unique = jadwiga_variants - dominik_variants
        
        total_unique_variants = len(dominik_variants | jadwiga_variants)
        jaccard_similarity = len(shared_variants) / total_unique_variants if total_unique_variants > 0 else 0
        
        comparison = {
            'dominik_variants': len(dominik_variants),
            'jadwiga_variants': len(jadwiga_variants),
            'shared_variants': len(shared_variants),
            'dominik_unique': len(dominik_unique),
            'jadwiga_unique': len(jadwiga_unique),
            'jaccard_similarity': jaccard_similarity,
            'overlap_percentage': (len(shared_variants) / min(len(dominik_variants), len(jadwiga_variants))) * 100 if min(len(dominik_variants), len(jadwiga_variants)) > 0 else 0
        }
        
        print(f"Dominik variants: {comparison['dominik_variants']:,}")
        print(f"Jadwiga variants: {comparison['jadwiga_variants']:,}")
        print(f"Shared variants: {comparison['shared_variants']:,}")
        print(f"Jaccard similarity: {comparison['jaccard_similarity']:.4f}")
        print(f"Overlap percentage: {comparison['overlap_percentage']:.2f}%")
        
        return comparison
    
    except Exception as e:
        print(f"Error in variant comparison: {e}")
        return {'error': str(e)}


def estimate_relationship(variant_comparison, cram_stats):
    """Estimate genetic relationship based on variant similarity"""
    
    print("\n=== Relationship Estimation ===")
    
    if 'error' in variant_comparison or not variant_comparison:
        print("Cannot estimate relationship - insufficient variant data")
        return {}
    
    jaccard = variant_comparison.get('jaccard_similarity', 0)
    overlap_pct = variant_comparison.get('overlap_percentage', 0)
    
    # Relationship estimation based on genetic similarity
    # Note: These thresholds are approximate and would need validation with known relationships
    if jaccard > 0.45 or overlap_pct > 90:
        relationship = "Identical twins or very close relatives"
        confidence = "High"
    elif jaccard > 0.25 or overlap_pct > 70:
        relationship = "Close family (parent-child, siblings)"
        confidence = "High"
    elif jaccard > 0.15 or overlap_pct > 50:
        relationship = "Related (grandparent-grandchild, uncle-nephew, etc.)"
        confidence = "Medium"
    elif jaccard > 0.08 or overlap_pct > 30:
        relationship = "Distant relatives"
        confidence = "Medium"
    else:
        relationship = "Unrelated or very distantly related"
        confidence = "Low"
    
    # Additional analysis using CRAM statistics if available
    coverage_analysis = {}
    if 'Dominik' in cram_stats and 'Jadwiga' in cram_stats:
        dom_stats = cram_stats['Dominik'].get('chromosome_stats', [])
        jad_stats = cram_stats['Jadwiga'].get('chromosome_stats', [])
        
        if dom_stats and jad_stats:
            # Compare coverage patterns
            dom_coverage = {stat['chromosome']: stat['mapped_reads'] for stat in dom_stats}
            jad_coverage = {stat['chromosome']: stat['mapped_reads'] for stat in jad_stats}
            
            common_chroms = set(dom_coverage.keys()) & set(jad_coverage.keys())
            if common_chroms:
                coverage_correlation = calculate_correlation(
                    [dom_coverage[c] for c in common_chroms],
                    [jad_coverage[c] for c in common_chroms]
                )
                coverage_analysis['chromosome_coverage_correlation'] = coverage_correlation
    
    relationship_analysis = {
        'estimated_relationship': relationship,
        'confidence': confidence,
        'jaccard_similarity': jaccard,
        'overlap_percentage': overlap_pct,
        'coverage_analysis': coverage_analysis,
        'interpretation': f"Based on {variant_comparison['shared_variants']:,} shared variants out of {variant_comparison['dominik_variants']:,} and {variant_comparison['jadwiga_variants']:,} total variants respectively."
    }
    
    print(f"Estimated relationship: {relationship}")
    print(f"Confidence: {confidence}")
    print(f"Interpretation: {relationship_analysis['interpretation']}")
    
    return relationship_analysis


def calculate_correlation(x, y):
    """Calculate Pearson correlation coefficient"""
    if len(x) != len(y) or len(x) < 2:
        return 0
    
    import numpy as np
    return np.corrcoef(x, y)[0, 1] if not np.isnan(np.corrcoef(x, y)[0, 1]) else 0


def load_cram_stats():
    """Load existing CRAM statistics"""
    stats_file = "/Users/dom/work/genetics/all_cram_stats.json"
    
    if os.path.exists(stats_file):
        try:
            with open(stats_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading CRAM stats: {e}")
    
    return {}


def main():
    """Main analysis function"""
    
    print("=== Simple CRAM and Genetic Analysis ===")
    print("Analyzing Dominik and Jadwiga's genetic data")
    
    # 1. Analyze existing VCF data
    print("\n" + "="*50)
    print("STEP 1: Analyzing VCF Data")
    print("="*50)
    
    vcf_analysis = analyze_existing_vcf_data()
    
    # 2. Compare variants
    print("\n" + "="*50)
    print("STEP 2: Variant Comparison")
    print("="*50)
    
    dominik_data = vcf_analysis.get('Dominik', {})
    jadwiga_data = vcf_analysis.get('Jadwiga', {})
    
    variant_comparison = compare_variants(dominik_data, jadwiga_data)
    
    # 3. Load CRAM statistics
    print("\n" + "="*50)
    print("STEP 3: CRAM Statistics")
    print("="*50)
    
    cram_stats = load_cram_stats()
    if cram_stats:
        print("Loaded existing CRAM statistics")
        for sample, stats in cram_stats.items():
            chrom_stats = stats.get('chromosome_stats', [])
            if chrom_stats:
                total_reads = sum(stat['mapped_reads'] for stat in chrom_stats)
                print(f"  {sample}: {total_reads:,} total mapped reads across {len(chrom_stats)} chromosomes/contigs")
    else:
        print("No CRAM statistics available")
    
    # 4. Estimate relationship
    print("\n" + "="*50)
    print("STEP 4: Relationship Analysis")
    print("="*50)
    
    relationship = estimate_relationship(variant_comparison, cram_stats)
    
    # 5. Compile and save results
    final_results = {
        'analysis_date': pd.Timestamp.now().isoformat(),
        'samples': ['Dominik', 'Jadwiga'],
        'vcf_analysis': vcf_analysis,
        'variant_comparison': variant_comparison,
        'cram_statistics': cram_stats,
        'relationship_analysis': relationship
    }
    
    # Save results
    output_file = "genetic_analysis_results.json"
    with open(output_file, 'w') as f:
        json.dump(final_results, f, indent=2, default=str)
    
    print(f"\n" + "="*50)
    print("ANALYSIS COMPLETE")
    print("="*50)
    print(f"Detailed results saved to: {output_file}")
    
    # Print summary
    if relationship:
        print(f"\nSUMMARY:")
        print(f"Estimated relationship: {relationship.get('estimated_relationship', 'Unknown')}")
        print(f"Confidence: {relationship.get('confidence', 'Unknown')}")
        if 'shared_variants' in variant_comparison:
            print(f"Shared variants: {variant_comparison['shared_variants']:,}")
            print(f"Genetic similarity: {variant_comparison['jaccard_similarity']:.4f}")


if __name__ == "__main__":
    main()



