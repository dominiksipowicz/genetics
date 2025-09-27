import itertools
from cyvcf2 import VCF

# Paths to your VCFs
dom_vcf_file = 'Dom/NG173LPFBH.vcf'
anna_vcf_file = 'Anna/NG1ABXTVKT.hard-filtered.vcf'
babcia_vcf_file = 'Babcia/M6FQWVFCC.hard-filtered.vcf'
jadwiga_vcf_file = 'Jadwiga/MG46K47NL.hard-filtered.vcf'


def variant_key(v):
    """Return CHROM:POS:REF>ALT for biallelic SNPs that PASS filters."""
    if v.FILTER is not None:  # skip non-PASS
        return None
    if len(v.ALT) != 1:
        return None
    ref, alt = v.REF, v.ALT[0]
    if len(ref) == 1 and len(alt) == 1:  # SNP only
        return f"{v.CHROM}:{v.POS}:{ref}>{alt}"
    return None


def extract_keys(vcf_file):
    keys = set()
    vcf = VCF(vcf_file)
    for v in vcf:
        k = variant_key(v)
        if k:
            keys.add(k)
    return keys


# Load sets
dom_keys = extract_keys(dom_vcf_file)
anna_keys = extract_keys(anna_vcf_file)
babcia_keys = extract_keys(babcia_vcf_file)
jadwiga_keys = extract_keys(jadwiga_vcf_file)

samples = {
    "Dom": dom_keys,
    "Anna": anna_keys,
    "Babcia": babcia_keys,
    "Jadwiga": jadwiga_keys
}

# Report totals
print("=== Total PASS biallelic SNPs ===")
for name, keys in samples.items():
    print(f"{name}: {len(keys):,}")

print("\n=== Pairwise overlaps (counts & Jaccard %) ===")
for (a, b) in itertools.combinations(samples.keys(), 2):
    set_a, set_b = samples[a], samples[b]
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    jaccard = inter / union * 100 if union else 0
    pct_a = inter / len(set_a) * 100 if set_a else 0
    pct_b = inter / len(set_b) * 100 if set_b else 0
    print(f"{a} vs {b}: shared {inter:,} | Jaccard {jaccard:.2f}% "
          f"| {pct_a:.2f}% of {a}, {pct_b:.2f}% of {b}")
