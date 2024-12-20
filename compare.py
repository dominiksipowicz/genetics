import pandas as pd
from cyvcf2 import VCF
from collections import Counter

# Paths to your and Anna's VCF files
dom_vcf_file = 'Dom/NG173LPFBH.vcf'
anna_vcf_file = 'Anna/NG1ABXTVKT.hard-filtered.vcf'
babcia_vcf_file = 'Babcia/M6FQWVFCC.hard-filtered.vcf'


def extract_variant_stats(vcf_file):
    vcf = VCF(vcf_file)
    ids = set()
    filter_stats = Counter()
    variant_stats = Counter({'SNP': 0, 'Indel': 0, 'SV': 0})

    for variant in vcf:
        if variant.ID:
            filter_stats[variant.FILTER] += 1
            if variant.FILTER is None:
                unique_id = f"{variant.CHROM}_{variant.POS}_{variant.ID}"
                ids.add(unique_id)
                ref_len = len(variant.REF)
                alt_len = len(variant.ALT[0])

                if 'SVTYPE' in variant.INFO:
                    variant_stats['SV'] += 1
                elif ref_len == 1 and alt_len == 1:
                    variant_stats['SNP'] += 1
                else:
                    variant_stats['Indel'] += 1
    
    return ids, filter_stats, variant_stats


# Extract variant IDs, filter statistics, and variant types from the VCF files
dom_ids, dom_filter_stats, dom_variant_stats = extract_variant_stats(dom_vcf_file)
anna_ids, anna_filter_stats, anna_variant_stats = extract_variant_stats(anna_vcf_file)
babcia_ids, babcia_filter_stats, babcia_variant_stats = extract_variant_stats(babcia_vcf_file)

# Find shared and unique IDs
dom_anna_shared_ids = dom_ids.intersection(anna_ids)
dom_babcia_shared_ids = dom_ids.intersection(babcia_ids)
anna_babcia_shared_ids = anna_ids.intersection(babcia_ids)

dom_anna_shared_ids_percent = len(dom_anna_shared_ids) / len(dom_ids) * 100
dom_babcia_shared_ids_percent = len(dom_babcia_shared_ids) / len(dom_ids) * 100
anna_babcia_shared_ids_percent = len(anna_babcia_shared_ids) / len(anna_ids) * 100

# Calculate the percentage of shared IDs based on Anna's total IDs
total_dom_ids = len(dom_ids)
total_anna_ids = len(anna_ids)
total_babcia_ids = len(babcia_ids)

# Save results to text files
with open('Dom/anna_dom_shared_ids.txt', 'w') as f:
    for id in dom_anna_shared_ids:
        f.write(id + '\n')

with open('Dom/babcia_dom_shared_ids.txt', 'w') as f:
    for id in dom_babcia_shared_ids:
        f.write(id + '\n')

with open('Anna/babcia_anna_shared_ids.txt', 'w') as f:
    for id in anna_babcia_shared_ids:
        f.write(id + '\n')

# Print variant statistics
print("\nVariant Statistics for Dom:")
for variant_type, count in dom_variant_stats.items():
    print(f"{variant_type}: {count}")

print("\nVariant Statistics for Anna:")
for variant_type, count in anna_variant_stats.items():
    print(f"{variant_type}: {count}")

print("\nVariant Statistics for Babcia:")
for variant_type, count in babcia_variant_stats.items():
    print(f"{variant_type}: {count}")

# Print summary statistics
print(f"\n")
print(f"Dom total variants: {total_dom_ids}")
print(f"Anna's total variants: {total_anna_ids}")
print(f"Babcia's total variants: {total_babcia_ids}")
print(f"Dom/Anna Shared IDs: {len(dom_anna_shared_ids)}")
print(f"Dom/Anna Percentage of shared IDs: {dom_anna_shared_ids_percent:.2f}%")
print(f"Dom/Babcia Shared IDs: {len(dom_babcia_shared_ids)}")
print(f"Dom/Babcia Percentage of shared IDs: {dom_babcia_shared_ids_percent:.2f}%")
print(f"Anna/Babcia Shared IDs: {len(anna_babcia_shared_ids)}")
print(f"Anna/Babcia Percentage of shared IDs: {anna_babcia_shared_ids_percent:.2f}%")
