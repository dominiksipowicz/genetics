#!/usr/bin/env python3
"""
Family Genetics Analysis - Proper Relatedness Calculation

This script performs proper genetic relatedness analysis by:
1. Normalizing VCFs using bcftools
2. Merging into multi-sample VCF  
3. Using PLINK2 for kinship coefficient calculation
4. Converting to % genome shared

Expected relationships:
- Jadwiga (mother) → Dom: ~50% shared
- Dom → Anna (daughter): ~50% shared  
- Jadwiga → Anna (grandmother): ~25% shared
- Babcia → Anna (great-grandmother): ~12.5% shared
"""

import subprocess
import os
import pandas as pd
from pathlib import Path

# Configuration
REFERENCE_FA = "Dom/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
WORK_DIR = "genetics_analysis"
PLINK2_PATH = "./plink2"  # Local PLINK2 binary

# VCF file paths (using compressed versions for better performance)
VCF_FILES = {
    "Dom": "Dom/NG173LPFBH.vcf.gz",
    "Anna": "Anna/NG1ABXTVKT.hard-filtered.vcf.gz", 
    "Babcia": "Babcia/M6FQWVFCC.hard-filtered.vcf.gz",
    "Jadwiga": "Jadwiga/MG46K47NL.hard-filtered.vcf.gz"
}

def run_cmd(cmd, description=""):
    """Run shell command with error handling."""
    print(f"Running: {description}")
    print(f"Command: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"ERROR in {description}:")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        raise RuntimeError(f"Command failed: {cmd}")
    
    print(f"✓ Completed: {description}\n")
    return result

def check_dependencies():
    """Check if required tools are installed."""
    missing = []
    
    # Check bcftools
    try:
        subprocess.run(["bcftools", "--version"], capture_output=True, check=True)
        print("✓ bcftools is installed")
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing.append("bcftools")
        print("✗ bcftools is missing")
    
    # Check local PLINK2 binary
    try:
        subprocess.run([PLINK2_PATH, "--version"], capture_output=True, check=True)
        print("✓ plink2 is available (local binary)")
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing.append("plink2")
        print("✗ plink2 local binary is missing or not executable")
    
    if missing:
        print("\nInstall missing tools:")
        if "bcftools" in missing:
            print("  brew install bcftools")
        if "plink2" in missing:
            print("  Local plink2 binary not found at:", PLINK2_PATH)
            print("  Make sure the plink2 binary is in the current directory")
        raise RuntimeError("Missing required tools")
    
    print("✓ All dependencies available\n")

def setup_work_directory():
    """Create and setup working directory."""
    Path(WORK_DIR).mkdir(exist_ok=True)
    print(f"✓ Work directory: {WORK_DIR}\n")

def normalize_vcfs():
    """Normalize all VCFs using bcftools norm."""
    print("=== STEP 1: Normalizing VCFs ===")
    
    normalized_files = {}
    
    for sample, vcf_path in VCF_FILES.items():
        if not os.path.exists(vcf_path):
            print(f"WARNING: {vcf_path} not found, skipping {sample}")
            continue
            
        # Output paths
        norm_vcf = f"{WORK_DIR}/{sample}_normalized.vcf.gz"
        
        # Create chromosome mapping file to remove "chr" prefix
        chr_map_file = f"{WORK_DIR}/chr_mapping.txt"
        if not os.path.exists(chr_map_file):
            with open(chr_map_file, 'w') as f:
                for i in range(1, 23):
                    f.write(f"chr{i}\t{i}\n")
                f.write("chrX\tX\n")
                f.write("chrY\tY\n") 
                f.write("chrM\tMT\n")
        
        # Process in steps to handle uncompressed VCF files properly
        temp_vcf = f"{WORK_DIR}/{sample}_temp.vcf.gz"
        
        # Step 1: Filter to main chromosomes and compress
        cmd = f"bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY {vcf_path} -Oz -o {temp_vcf}"
        run_cmd(cmd, f"Filtering {sample} to main chromosomes")
        
        # Step 2: Index the temp file
        run_cmd(f"bcftools index {temp_vcf}", f"Indexing temp {sample} VCF")
        
        # Step 3: Rename chromosomes and normalize
        cmd = f"""bcftools annotate --rename-chrs {chr_map_file} {temp_vcf} | \
                  bcftools norm -m-any -f {REFERENCE_FA} | \
                  bcftools norm -d all -f {REFERENCE_FA} -Oz -o {norm_vcf}"""
        
        run_cmd(cmd, f"Normalizing {sample} VCF")
        
        # Index the normalized VCF
        run_cmd(f"bcftools index {norm_vcf}", f"Indexing {sample} normalized VCF")
        
        # Clean up temp file and its index (could be .csi or .tbi)
        os.remove(temp_vcf)
        for ext in ['.csi', '.tbi']:
            temp_idx = f"{temp_vcf}{ext}"
            if os.path.exists(temp_idx):
                os.remove(temp_idx)
        
        normalized_files[sample] = norm_vcf
    
    return normalized_files

def merge_vcfs(normalized_files):
    """Merge normalized VCFs into multi-sample file."""
    print("=== STEP 2: Merging VCFs ===")
    
    vcf_list = " ".join(normalized_files.values())
    merged_vcf = f"{WORK_DIR}/family_merged.vcf.gz"
    
    # Merge VCFs
    cmd = f"bcftools merge {vcf_list} -Oz -o {merged_vcf}"
    run_cmd(cmd, "Merging family VCFs")
    
    # Index merged VCF
    run_cmd(f"bcftools index {merged_vcf}", "Indexing merged VCF")
    
    return merged_vcf

def filter_for_relatedness(merged_vcf):
    """Filter VCF for high-quality sites suitable for relatedness analysis."""
    print("=== STEP 3: Filtering for Relatedness Analysis ===")
    
    filtered_vcf = f"{WORK_DIR}/family_filtered.vcf.gz"
    
    # Filter to PASS, biallelic SNPs, autosomes only, with good genotypes
    filter_expr = """FILTER="PASS" && TYPE="snp" && N_ALT=1 && 
                     CHROM!="X" && CHROM!="Y" && CHROM!="MT" &&
                     F_MISSING < 0.1"""
    
    cmd = f"""bcftools view -i '{filter_expr}' {merged_vcf} -Oz -o {filtered_vcf}"""
    run_cmd(cmd, "Filtering VCF for relatedness analysis")
    
    # Index filtered VCF  
    run_cmd(f"bcftools index {filtered_vcf}", "Indexing filtered VCF")
    
    return filtered_vcf

def calculate_kinship(filtered_vcf):
    """Use PLINK2 to calculate kinship coefficients."""
    print("=== STEP 4: Calculating Kinship Coefficients ===")
    
    # Convert VCF to PLINK format
    plink_prefix = f"{WORK_DIR}/family"
    
    cmd = f"{PLINK2_PATH} --vcf {filtered_vcf} --make-bed --out {plink_prefix}"
    run_cmd(cmd, "Converting VCF to PLINK format")
    
    # Calculate kinship matrix
    cmd = f"{PLINK2_PATH} --bfile {plink_prefix} --make-king-table --out {plink_prefix}"
    run_cmd(cmd, "Calculating kinship coefficients")
    
    return f"{plink_prefix}.kin0"

def analyze_kinship_results(king_file):
    """Parse PLINK2 kinship results and convert to % genome shared."""
    print("=== STEP 5: Analyzing Kinship Results ===")
    
    # Read PLINK2 kinship table
    df = pd.read_csv(king_file, sep='\t')
    
    print("Raw kinship coefficients:")
    print(df[['IID1', 'IID2', 'KINSHIP']].to_string(index=False))
    
    print("\n=== Family Relatedness Analysis ===")
    print("Expected vs Observed Kinship Coefficients:")
    print("(Kinship coefficient ≈ % genome shared / 100)\n")
    
    # Expected relationships (using actual sample IDs from VCF files)
    expected = {
        ('NG173LPFBH', 'NG1ABXTVKT'): {'kinship': 0.25, 'relationship': 'Parent-Child (Dom-Anna)', 'expected_pct': 50},
        ('MG46K47NL', 'NG173LPFBH'): {'kinship': 0.25, 'relationship': 'Parent-Child (Jadwiga-Dom)', 'expected_pct': 50}, 
        ('MG46K47NL', 'NG1ABXTVKT'): {'kinship': 0.125, 'relationship': 'Grandparent-Grandchild (Jadwiga-Anna)', 'expected_pct': 25},
        ('M6FQWVFCC', 'NG1ABXTVKT'): {'kinship': 0.0625, 'relationship': 'Great-Grandparent-Great-Grandchild (Babcia-Anna)', 'expected_pct': 12.5},
        ('M6FQWVFCC', 'NG173LPFBH'): {'kinship': 0.125, 'relationship': 'Grandparent-Grandchild (Babcia-Dom)', 'expected_pct': 25},
        ('M6FQWVFCC', 'MG46K47NL'): {'kinship': 0.25, 'relationship': 'Parent-Child (Babcia-Jadwiga)', 'expected_pct': 50}
    }
    
    results = []
    
    for _, row in df.iterrows():
        id1, id2 = sorted([row['IID1'], row['IID2']])
        observed_kinship = row['KINSHIP']
        observed_pct = observed_kinship * 100 * 2  # Convert kinship to % genome shared
        
        pair = (id1, id2)
        if pair in expected:
            exp = expected[pair]
            diff = abs(observed_pct - exp['expected_pct'])
            
            print(f"{id1} ↔ {id2}:")
            print(f"  Relationship: {exp['relationship']}")
            print(f"  Expected: {exp['expected_pct']:.1f}% genome shared")
            print(f"  Observed: {observed_pct:.1f}% genome shared") 
            print(f"  Difference: {diff:.1f}%")
            print(f"  Kinship coefficient: {observed_kinship:.4f}")
            print()
            
            results.append({
                'Pair': f"{id1}-{id2}",
                'Relationship': exp['relationship'],
                'Expected_Pct': exp['expected_pct'],
                'Observed_Pct': round(observed_pct, 1),
                'Difference': round(diff, 1),
                'Kinship': round(observed_kinship, 4)
            })
    
    # Save results to CSV
    results_df = pd.DataFrame(results)
    results_file = f"{WORK_DIR}/family_relatedness_results.csv"
    results_df.to_csv(results_file, index=False)
    print(f"✓ Results saved to: {results_file}")
    
    # Answer the main question
    print("\n" + "="*60)
    print("ANSWER TO YOUR QUESTION:")
    print("How much % of you (Dom) did Anna inherit?")
    
    dom_anna = results_df[results_df['Pair'] == 'NG173LPFBH-NG1ABXTVKT']
    if not dom_anna.empty:
        pct = dom_anna['Observed_Pct'].iloc[0]
        print(f"Anna inherited approximately {pct}% of your genome.")
        print(f"This is {'normal' if abs(pct - 50) < 5 else 'unusual'} for a parent-child relationship.")
    
    print("="*60)
    
    return results_df

def main():
    """Main analysis pipeline."""
    print("Family Genetics Analysis - Proper Relatedness Calculation")
    print("="*60)
    
    try:
        # Setup
        check_dependencies()
        setup_work_directory()
        
        # Analysis pipeline
        normalized_files = normalize_vcfs()
        merged_vcf = merge_vcfs(normalized_files)
        filtered_vcf = filter_for_relatedness(merged_vcf)
        king_file = calculate_kinship(filtered_vcf)
        results = analyze_kinship_results(king_file)
        
        print("\n✓ Analysis completed successfully!")
        print(f"All intermediate files saved in: {WORK_DIR}/")
        
    except Exception as e:
        print(f"\n✗ Analysis failed: {e}")
        raise

if __name__ == "__main__":
    main()
