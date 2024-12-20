import pandas as pd
from cyvcf2 import VCF
from collections import Counter

# Path to your VCF file
# Dom
# vcf_file = 'Dom/NG173LPFBH.vcf'
# output_file = 'Dom/parsed_vcf_data.csv'
# Anna
vcf_file = 'Anna/NG1ABXTVKT.hard-filtered.vcf'
output_file = 'Anna/parsed_vcf_data.csv'
# Babcia
# vcf_file = 'Babcia/M6FQWVFCC.hard-filtered.vcf'
# output_file = 'Babcia/parsed_vcf_data.csv'

# Function to parse VCF and extract relevant data
def parse_vcf(vcf_file):
    vcf = VCF(vcf_file)
    records = []
    
    for variant in vcf:
        for alt in variant.ALT:
            gene_info = {
                'CHROM': variant.CHROM,
                'POS': variant.POS,
                'ID': variant.ID,
                'REF': variant.REF,
                'ALT': alt,
                'QUAL': variant.QUAL,
                'FILTER': variant.FILTER,
                'INFO': variant.INFO,
                'FORMAT': variant.FORMAT,
                'SAMPLE': variant.genotypes[0]
            }
            records.append(gene_info)
    
    return pd.DataFrame(records)


# Parse the VCF file
try:
    df = parse_vcf(vcf_file)
except Exception as e:
    print(f"Error reading VCF file: {e}")
    exit(1)

# Print the dataframe
print("\nGene Data:")
print(df.head())  # Print the first 5 rows for brevity

# Perform a simple analysis: count the number of entries per chromosome
chromosome_counts = df['CHROM'].value_counts()

# Print the analysis result
print("\nNumber of entries per chromosome:")
print(chromosome_counts)

# Additional analysis: count the number of unique mutations (ALT) and the most frequent mutation
unique_mutations = df['ALT'].nunique()
most_frequent_mutation = df['ALT'].value_counts().idxmax()

# Average QUAL
average_qual = df['QUAL'].mean()

# Filter value counts
filter_counts = df['FILTER'].value_counts()

print(f"\nNumber of unique mutations: {unique_mutations}")
print(f"Most frequent mutation: {most_frequent_mutation}")
print(f"Average QUAL: {average_qual}")
print(f"Filter value counts: {filter_counts}")

# Save the dataframe to a CSV file for further analysis if needed
df.to_csv(output_file, index=False)
