import os
import pandas as pd
import glob
import sys

# ================= CONFIGURATION =================
# Root directory where script 02 saves the reports
REPORTS_ROOT = "../results/reports/"
OUTPUT_DIR = "../results/final_tables/"

# Ranks to extract
TARGET_RANKS = {
    'P': 'Phylum',
    'C': 'Class',
    'O': 'Order',
    'F': 'Family',
    'G': 'Genus',
    'S': 'Species'
}
# =================================================

def get_db_selection():
    """
    Scans the report directory and asks the user which Database to process.
    """
    # 1. Automatic: Check command line arg (e.g., python script.py FISH)
    if len(sys.argv) > 1:
        return sys.argv[1]

    # 2. Interactive: Scan folders
    if not os.path.exists(REPORTS_ROOT):
        print(f"❌ Error: Directory {REPORTS_ROOT} does not exist.")
        sys.exit(1)

    options = [d for d in os.listdir(REPORTS_ROOT) 
               if os.path.isdir(os.path.join(REPORTS_ROOT, d))]
    
    if not options:
        print(f"❌ No report folders found in {REPORTS_ROOT}")
        sys.exit(1)

    print("\n📊 AVAILABLE DATABASES:")
    for i, db in enumerate(options):
        print(f"   [{i + 1}] {db}")

    while True:
        try:
            choice = input("\n👉 Select the database number: ")
            choice_idx = int(choice) - 1
            if 0 <= choice_idx < len(options):
                return options[choice_idx]
            print("⚠️  Invalid number. Try again.")
        except ValueError:
            print("⚠️  Please enter a number.")

def get_metric_selection():
    """
    Asks the user which data column to extract from the Kraken report.
    Returns: (column_index, metric_label)
    """
    print("\n📉 DATA METRIC SELECTION:")
    print("   [1] Cumulative/Rooted Reads (Includes reads from sub-species)")
    print("   [2] Direct Reads (Reads assigned specifically to this rank)")
    
    # Check if a second argument was passed in command line for automation
    # usage: python script.py DB_NAME METRIC_ID
    if len(sys.argv) > 2:
        choice = sys.argv[2]
        print(f"🔹 Automatic Selection: Option {choice}")
    else:
        choice = input("\n👉 Select the metric to extract (1 or 2): ")

    if choice == '2':
        # parts[2] is Direct Reads
        return 2, "Direct_Reads"
    else:
        # Default to parts[1] (Cumulative) if 1 or anything else is selected
        # This matches your original "cite 1" logic
        return 1, "Cumulative_Reads"

def parse_kraken_report(filepath, sample_name, col_index):
    """
    Parses a single report extracting the specific column index selected by user.
    """
    extracted_data = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6: continue
                
                try:
                    # MODULARITY: Uses the col_index passed from user selection
                    # parts[1] = Rooted, parts[2] = Direct
                    count = int(parts[col_index])  
                    
                    rank_code = parts[3]
                    tax_id = int(parts[4])
                    name = parts[5].strip()
                except ValueError:
                    continue

                if rank_code in TARGET_RANKS:
                    extracted_data.append({
                        'Rank': TARGET_RANKS[rank_code],
                        'TaxID': tax_id,
                        'Name': name,
                        'Sample': sample_name,
                        'Reads': count
                    })
    except Exception as e:
        print(f"❌ Error reading {filepath}: {e}")
        
    return extracted_data

def main():
    print("==========================================")
    print("🚀 KRAKEN2 TABLE GENERATOR")
    print("==========================================")
    
    # --- PHASE 1: SELECTION ---
    selected_db = get_db_selection()
    col_index, metric_name = get_metric_selection()
    
    input_path = os.path.join(REPORTS_ROOT, selected_db)
    
    # Output file now includes BOTH the DB name and the Metric name
    # Example: Taxonomy_FISH_Cumulative_Reads.xlsx
    output_filename = f"Taxonomy_{selected_db}_{metric_name}.xlsx"
    output_file = os.path.join(OUTPUT_DIR, output_filename)

    print(f"\n📂 Processing Folder: {input_path}")
    print(f"🔍 Extraction Mode:   {metric_name.replace('_', ' ')}")
    
    # --- PHASE 2: PROCESSING ---
    report_files = glob.glob(os.path.join(input_path, "*_report.txt"))
    
    if not report_files:
        print(f"❌ No report files found in {input_path}")
        return

    print(f"   └── Found {len(report_files)} samples.")

    all_data = []
    for report in report_files:
        filename = os.path.basename(report)
        sample_name = filename.replace('_report.txt', '')
        
        # Pass the col_index to the parser
        data = parse_kraken_report(report, sample_name, col_index)
        all_data.extend(data)

    if not all_data:
        print("⚠️  No data found matching criteria.")
        return

    # --- PHASE 3: SAVING ---
    print("🔄 Merging and pivoting data...")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    df = pd.DataFrame(all_data)
    
    # Pivot Table
    pivot_df = df.pivot_table(
        index=['Rank', 'TaxID', 'Name'], 
        columns='Sample', 
        values='Reads', 
        fill_value=0
    ).reset_index()

    # Sorting
    rank_order = list(TARGET_RANKS.values())
    pivot_df['Rank_Cat'] = pd.Categorical(pivot_df['Rank'], categories=rank_order, ordered=True)
    pivot_df = pivot_df.sort_values(['Rank_Cat', 'Name']).drop(columns=['Rank_Cat'])

    try:
        pivot_df.to_excel(output_file, index=False)
        print("\n🎉 SUCCESS!")
        print(f"💾 Table saved to: {output_file}")
    except Exception as e:
        print(f"❌ Error saving Excel: {e}")

if __name__ == "__main__":
    main()