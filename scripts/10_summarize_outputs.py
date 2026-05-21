import pandas as pd
import os
import glob
import argparse

def process_excel_file(file_path):
    """
    Processes a single Excel file.
    Returns a DataFrame with Unclassified, Classified, and Total reads for each sample.
    """
    print(f"Processing: {file_path}")
    
    # Read the excel file
    df = pd.read_excel(file_path)
    
    # Identify sample columns dynamically (assuming first 3 are Rank, TaxID, Name)
    sample_cols = df.columns[3:] 
    
    summary_data = []
    
    for col in sample_cols:
        # Isolate the unclassified row robustly (lowercased to avoid case sensitivity issues)
        unclass_rows = df[df['Rank'].astype(str).str.lower() == 'unclassified']
        
        if not unclass_rows.empty:
            unclassified = unclass_rows[col].values[0]
        else:
            unclassified = 0
            
        total = df[col].sum()
        classified = total - unclassified
        
        summary_data.append({
            'Source_DB': os.path.basename(file_path),
            'Sample': col,
            'Unclassified': unclassified,
            'Classified': classified,
            'Total': total
        })
        
    return pd.DataFrame(summary_data)


def main(input_dir, output_dir):
    """
    Iterates through all .xlsx files, creates individual summaries,
    and calculates the grand total matrix.
    """
    # 1. Setup Output Directories
    os.makedirs(output_dir, exist_ok=True)
    per_file_dir = os.path.join(output_dir, 'per_db_summaries')
    os.makedirs(per_file_dir, exist_ok=True)
    
    # 2. Locate all Excel files
    search_pattern = os.path.join(input_dir, "*.xlsx")
    excel_files = glob.glob(search_pattern)
    
    if not excel_files:
        print(f"Error: No .xlsx files found in '{input_dir}'")
        return

    all_summaries = []

    # 3. Process each file modularly
    for file in excel_files:
        # Ignore temporary Excel locks
        if os.path.basename(file).startswith('~$'):
            continue
            
        # Process the single dataframe
        df_summary = process_excel_file(file)
        all_summaries.append(df_summary)
        
        # Save the independent DB summary
        base_name = os.path.splitext(os.path.basename(file))[0]
        out_file = os.path.join(per_file_dir, f"{base_name}_summary.csv")
        df_summary.to_csv(out_file, index=False)
        print(f"  --> Saved individual DB summary: {out_file}")

    # 4. Aggregate Master Matrix
    combined_df = pd.concat(all_summaries, ignore_index=True)
    
    # Calculate the Grand Total across ALL samples in ALL DBs
    grand_total = combined_df[['Unclassified', 'Classified', 'Total']].sum()
    grand_total['Source_DB'] = 'ALL_DATABASES'
    grand_total['Sample'] = 'GRAND_TOTAL'
    
    # Append the Grand Total vector to the matrix
    combined_df.loc[len(combined_df)] = grand_total
    
    # 5. Export Master Table
    final_out_path = os.path.join(output_dir, "master_classification_total.csv")
    combined_df.to_csv(final_out_path, index=False)
    print(f"\nSuccess! Master matrix saved to: {final_out_path}")


if __name__ == "__main__":
    # Configure your paths here or run via CLI
    parser = argparse.ArgumentParser(description="Kraken/Taxa Read Summarizer")
    parser.add_argument("-i", "--input", default=".", help="Directory containing .xlsx files")
    parser.add_argument("-o", "--output", default="result/final_tables", help="Output directory")
    
    args = parser.parse_args()
    
    main(args.input, args.output)
