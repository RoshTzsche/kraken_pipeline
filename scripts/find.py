import pandas as pd

# Load your metadata IDs
meta_df = pd.read_excel("../data/metadata.xlsx")
meta_samples = set(meta_df['SampleID'].astype(str).str.strip().str.lower())

# Load your data table columns
data_df = pd.read_excel("../results/final_tables/Taxonomy_PLANTS_Cumulative_Reads.xlsx", sheet_name=0)
# Filter out the annotation columns to just get sample names
data_samples = set([col.split('_')[0].strip().lower() for col in data_df.columns 
                    if col not in ['Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name']])

# Find the missing sample
missing = meta_samples - data_samples
print(f"Sample in metadata but missing from data: {missing}")