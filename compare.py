import pandas as pd
from cyvcf2 import VCF
from collections import Counter

def extract_passed_ids_and_filter_stats(vcf_file):
    vcf = VCF(vcf_file)
    ids = set()
    filter_stats = Counter()

    for variant in vcf:
        if variant.ID:
            filter_stats[variant.FILTER] += 1
            if variant.FILTER is None:
                unique_id = f"{variant.CHROM}_{variant.POS}_{variant.ID}"
                ids.add(unique_id)
    
    return ids, filter_stats

# Paths to your and Anna's VCF files
my_vcf_file = 'Dom/NG173LPFBH.vcf'
anna_vcf_file = 'Anna/NG1ABXTVKT.hard-filtered.vcf'

# Extract variant IDs and filter statistics from the VCF files
my_ids, my_filter_stats = extract_passed_ids_and_filter_stats(my_vcf_file)
anna_ids, anna_filter_stats = extract_passed_ids_and_filter_stats(anna_vcf_file)

# Find shared and unique IDs
shared_ids = my_ids.intersection(anna_ids)
my_unique_ids = my_ids.difference(anna_ids)
anna_unique_ids = anna_ids.difference(my_ids)

# Calculate the percentage of shared IDs based on Anna's total IDs
total_my_ids = len(my_ids)
total_anna_ids = len(anna_ids)
shared_ids_percent = len(shared_ids) / total_anna_ids * 100

# Save results to text files
with open('Dom/shared_ids.txt', 'w') as f:
    for id in shared_ids:
        f.write(id + '\n')

with open('Dom/my_unique_ids.txt', 'w') as f:
    for id in my_unique_ids:
        f.write(id + '\n')

with open('Dom/anna_unique_ids.txt', 'w') as f:
    for id in anna_unique_ids:
        f.write(id + '\n')

# Print filter statistics
print("FILTER Statistics for Myself:")
for filter_value, count in my_filter_stats.items():
    print(f"{filter_value}: {count}")

print("\nFILTER Statistics for Anna:")
for filter_value, count in anna_filter_stats.items():
    print(f"{filter_value}: {count}")

# Print summary statistics
print(f"\nMy total variants: {total_my_ids}")
print(f"Anna's total variants: {total_anna_ids}")
print(f"Shared IDs: {len(shared_ids)}")
print(f"Percentage of shared IDs: {shared_ids_percent:.2f}%")
