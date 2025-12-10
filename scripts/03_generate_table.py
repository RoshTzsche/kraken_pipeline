import os
import pandas as pd
import glob

# ================= CONFIGURATION =================
# Paths relative to the scripts/ directory
INPUT_DIR = "../results/reports/"
OUTPUT_FILE = "../results/final_tables/Taxonomy_Stacked_Results.xlsx"

# Ranks to extract and their label in the final Excel
TARGET_RANKS = {
    'P': 'Phylum',
    'C': 'Class',
    'O': 'Order',
    'F': 'Family',
    'G': 'Genus',
    'S': 'Species'
}
# =================================================

def parse_kraken_report_for_ranks(filepath, sample_name):
    """
    Parses a Kraken2 report to extract TaxID, Name, and cumulative Reads
    for specific taxonomic ranks.
    """
    extracted_data = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6: continue
                
                try:
                    # Kraken format: % | reads_rooted | reads_direct | rank | taxid | name
                    [cite_start]count = int(parts[1]) # We use reads_rooted (cumulative) [cite: 1]
                    rank_code = parts[3]
                    tax_id = int(parts[4])
                    name = parts[5].strip()
                except ValueError:
                    continue

                if rank_code in TARGET_RANKS:
                    rank_name = TARGET_RANKS[rank_code]
                    
                    extracted_data.append({
                        'Rank': rank_name,
                        'TaxID': tax_id,
                        'Name': name,
                        'Sample': sample_name,
                        'Reads': count
                    })
    except Exception as e:
        print(f"❌ Error reading {filepath}: {e}")
        
    return extracted_data

def main():
    print("🚀 Generating consolidated taxonomy table...")
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    # 1. Find report files
    report_files = glob.glob(os.path.join(INPUT_DIR, "*_report.txt"))
    
    if not report_files:
        print(f"❌ No reports found in {INPUT_DIR}")
        print("   Did you run script '02_run_kraken.sh'?")
        return

    print(f"📂 Processing {len(report_files)} report files...")

    all_data = []

    # 2. Extract data from each report
    for report in report_files:
        filename = os.path.basename(report)
        # Simple cleanup (Sample_report.txt -> Sample)
        sample_name = filename.replace('_report.txt', '')
        
        print(f"  └── Reading: {sample_name}")
        data = parse_kraken_report_for_ranks(report, sample_name)
        all_data.extend(data)

    if not all_data:
        print("⚠️  No relevant taxonomic data found.")
        return

    # 3. Create DataFrame and Pivot
    print("🔄 Merging and structuring data...")
    df = pd.DataFrame(all_data)

    # Pivot: Index (Taxonomy) x Columns (Samples)
    pivot_df = df.pivot_table(
        index=['Rank', 'TaxID', 'Name'], 
        columns='Sample', 
        values='Reads', 
        fill_value=0
    ).reset_index()

    # 4. Hierarchical Sorting
    # Define categorical order for ranks to appear logically in the Excel
    rank_order = list(TARGET_RANKS.values()) # ['Phylum', 'Class', ...]
    pivot_df['Rank_Cat'] = pd.Categorical(pivot_df['Rank'], categories=rank_order, ordered=True)
    
    # Sort by Rank Category and then by Name
    pivot_df = pivot_df.sort_values(['Rank_Cat', 'Name']).drop(columns=['Rank_Cat'])

    # 5. Save to Excel
    try:
        pivot_df.to_excel(OUTPUT_FILE, index=False)
        print(f"💾 Success! Table saved to:\n   {OUTPUT_FILE}")
    except Exception as e:
        print(f"❌ Error saving Excel file: {e}")

if __name__ == "__main__":
    main()