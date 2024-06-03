import pandas as pd

# Load the gene data from a CSV file with error handling
try:
    df = pd.read_csv('Dom/NG173LPFBH.csv', on_bad_lines='warn')
except pd.errors.ParserError as e:
    print(f"Error reading CSV file: {e}")
    exit(1)

# Print the dataframe
print("\nGene Data:")
print(df.head())  # Print the first 5 rows for brevity

# # Perform a simple analysis: count the number of entries per chromosome
# chromosome_counts = df['CHROM'].value_counts()

# # Print the analysis result
# print("\nNumber of entries per chromosome:")
# print(chromosome_counts)

# # Additional analysis: count the number of unique mutations (ALT) and the most frequent mutation
# unique_mutations = df['ALT'].nunique()
# most_frequent_mutation = df['ALT'].value_counts().idxmax()

# print(f"\nNumber of unique mutations: {unique_mutations}")
# print(f"Most frequent mutation: {most_frequent_mutation}")
