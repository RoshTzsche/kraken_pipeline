import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def generate_microbiome_pdf(data_path, rank_level='family', threshold=0.01, output_file='microbiome_report.pdf'):
    """
    Generates a PDF with a stacked bar chart of the relative abundance of the microbiome,
    including metadata of the parameters used directly in the chart.
    """
    
    # 1. Load and segment data by Taxonomic Level
    #df = pd.read_csv(data_path)
    df = pd.read_excel(data_path, sheet_name=0) 
    # 1. Load Data
    df = pd.read_excel(data_path, sheet_name=0)
    # Dynamically handle CSV or Excel metadata
    if metadata_path.lower().endswith('.csv'):
        meta_df = pd.read_csv(metadata_path)
    elif metadata_path.lower().endswith(('.xls', '.xlsx')):
        meta_df = pd.read_excel(metadata_path)
    else:
        raise ValueError("Metadata file must be a .csv or .xlsx")
 
    df_rank = df[df['Rank'].str.lower() == rank_level.lower()].copy()
    
    if df_rank.empty:
        raise ValueError(f"No data found for the taxonomic level: {rank_level}")
        
    df_rank.set_index('Scientific Name', inplace=True)
    
    # Isolate numeric columns from samples
    sample_cols = [col for col in df_rank.columns if col not in ['Rank', 'TaxId', 'original_header']]
    df_counts = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # 2. Transformation to Relative Abundance (p_{ij} = x_{ij} / sum(x_{*j}))
    rel_abund = df_counts.div(df_counts.sum(axis=0), axis=1)
    
    # 3. Mean calculation and Threshold Filtering (tau)
    mean_rel_abund = rel_abund.mean(axis=1)
    keep_taxa = mean_rel_abund[mean_rel_abund >= threshold].index
    
    df_kept = rel_abund.loc[keep_taxa]
    
    # Agglomeration of the "Others" category
    df_others_sum = rel_abund.loc[~rel_abund.index.isin(keep_taxa)].sum(axis=0)
    df_others = pd.DataFrame(df_others_sum, columns=['Others']).T
    
    # Concatenate and transpose for plotting (X-axis = Samples)
    df_final = pd.concat([df_kept, df_others])
    df_plot = df_final.T
    
    # Extended color palette
    colors = [
        "#b988d5", "#cbd588", "#88a05b", "#ffe156", "#6b8ba4", "#fda47d", 
        "#8e5634", "#d44645", "#5d9ca3", "#63b7af", "#dcd3ff", "#ff94cc", 
        "#ffa45b", "#806d40", "#2a363b", "#99b898", "#feceab", "#ff847c", 
        "#e84a5f", "#2a363b", "#56a5cc", "#c3e88d", "#ffcc5c", "#b09f59", 
        "#ff5e57", "#674d3c", "#4c4f69", "#8372a8", "#ff7c43", "#6a8a82"
    ]
    
    # 4. PDF Creation and Vector Rendering
    with PdfPages(output_file) as pdf:
        # Define canvas dimensions (Width x Height)
        fig, ax = plt.subplots(figsize=(14, 10))
        
        # Assign colors, forcing "Others" to always be gray (#cccccc)
        plot_colors = colors[:len(keep_taxa)] + ['#cccccc']
        
        # Generate stacked bars
        df_plot.plot(kind='bar', stacked=True, ax=ax, width=0.85, color=plot_colors)
        
        # -- METADATA IN THE CHART --
        # Main title indicating the level (Family, Genus, etc.)
        plt.title(f'Microbiome Composition - Grouping: {rank_level.capitalize()}', 
                  fontsize=16, pad=25)
        
        # Subtitle (suptitle) indicating the mathematical Threshold (tau) and taxa count
        plt.suptitle(f'Filtering threshold ($\\tau$): {threshold*100}% | Taxa shown: {len(keep_taxa)} + Others', 
                     fontsize=12, color='black', y=0.88)
        
        # Axis names
        plt.ylabel('Relative Abundance', fontsize=12, fontweight='bold')
        plt.xlabel('Sample Identifier', fontsize=12, fontweight='bold')
        
        # Rotate X-axis names to prevent collision
        plt.xticks(rotation=45, ha='right', fontsize=9)
        
        # Locate the legend outside the chart area
        plt.legend(title=rank_level.capitalize(), 
                   bbox_to_anchor=(1.01, 1), # Coordinates: Starts just outside on X (1.01), and top on Y (1)
                   loc='upper left', 
                   fontsize=9, 
                   borderaxespad=0.)
        
        # tight_layout compresses the chart to respect the space for the legend and titles
        plt.tight_layout()
        
        # Save the current figure to the PDF document
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()


def generate_grouped_microbiome_pdf(data_path, metadata_path, category_col, sample_id_col='SampleID', rank_level='family', threshold=0.01, output_file='grouped_report.pdf'):
    """
    Groups absolute counts by a metadata category before calculating relative abundance.
    """
    # 1. Load Data
    df = pd.read_excel(data_path, sheet_name=0)
    # Dynamically handle CSV or Excel metadata
    if metadata_path.lower().endswith('.csv'):
        meta_df = pd.read_csv(metadata_path)
    elif metadata_path.lower().endswith(('.xls', '.xlsx')):
        meta_df = pd.read_excel(metadata_path)
    else:
        raise ValueError("Metadata file must be a .csv or .xlsx")
    meta_df = pd.read_excel(metadata_path, sheet_name=0) # Using read_csv as per your uploaded file
    
    print(f"Columns found in metadata: {meta_df.columns.tolist()}")
    df_rank = df[df['Rank'].str.lower() == rank_level.lower()].copy()
    if df_rank.empty:
        raise ValueError(f"No data found for the taxonomic level: {rank_level}")
        
    df_rank.set_index('Scientific Name', inplace=True)
    
    # 2. Isolate numeric columns (Samples)
    sample_cols = [col for col in df_rank.columns if col not in ['Rank', 'TaxId', 'original_header']]
    df_counts = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # 3. Map samples to categories and sum counts
    # Force metadata sample IDs to be strings and strip any hidden spaces
    meta_df[sample_id_col] = meta_df[sample_id_col].astype(str).str.strip()
    
    # Create a dictionary mapping the sample IDs to the requested category
    category_map = dict(zip(meta_df[sample_id_col], meta_df[category_col]))
    
    # Force taxonomy column names to be clean strings too
    df_counts.columns = [str(col).strip() for col in df_counts.columns]
    print("\n--- DEBUGGING SAMPLE IDs ---")
    print(f"First 5 IDs in Metadata: {meta_df[sample_id_col].tolist()[:5]}")
    print(f"First 5 Columns in Taxonomy: {list(df_counts.columns)[:5]}")
    print("----------------------------\n")
    # ----------------------------
    df_counts.columns = [str(col).rsplit('_', 1)[0].strip() if '_' in str(col) else str(col).strip() for col in df_counts.columns]
    cols_to_keep = [col for col in df_counts.columns if col in category_map]
    if not cols_to_keep:
        raise ValueError("CRITICAL ERROR: None of the Sample IDs in your data matched the metadata.")
    
    df_counts = df_counts[cols_to_keep]
    
    # Rename sample columns to their respective categories (e.g., "Male", "Female")
    df_counts_mapped = df_counts.rename(columns=category_map)
    
    # Group identical column names and sum their counts across axis=1
    df_grouped_counts = df_counts_mapped.groupby(df_counts_mapped.columns, axis=1).sum()
    # 4. Transformation to Relative Abundance ($ p_{i,k} $)
    rel_abund = df_grouped_counts.div(df_grouped_counts.sum(axis=0), axis=1)
    
    # 5. Mean calculation, Thresholding (\tau), and Agglomeration
    mean_rel_abund = rel_abund.mean(axis=1)
    keep_taxa = mean_rel_abund[mean_rel_abund >= threshold].index
    
    df_kept = rel_abund.loc[keep_taxa]
    df_others = pd.DataFrame(rel_abund.loc[~rel_abund.index.isin(keep_taxa)].sum(axis=0), columns=['Others']).T
    
    df_plot = pd.concat([df_kept, df_others]).T
    
    # Extended color palette
    colors = [
        "#b988d5", "#cbd588", "#88a05b", "#ffe156", "#6b8ba4", "#fda47d", 
        "#8e5634", "#d44645", "#5d9ca3", "#63b7af", "#dcd3ff", "#ff94cc", 
        "#ffa45b", "#806d40", "#2a363b", "#99b898", "#feceab", "#ff847c", 
        "#e84a5f", "#2a363b", "#56a5cc", "#c3e88d", "#ffcc5c", "#b09f59", 
        "#ff5e57", "#674d3c", "#4c4f69", "#8372a8", "#ff7c43", "#6a8a82"
    ]
    
    # 6. Render
    with PdfPages(output_file) as pdf:
        fig, ax = plt.subplots(figsize=(14, 10))
        plot_colors = colors[:len(keep_taxa)] + ['#cccccc']
        
        df_plot.plot(kind='bar', stacked=True, ax=ax, width=0.85, color=plot_colors)
        
        plt.title(f'Microbiome Composition by {category_col} - Grouping: {rank_level.capitalize()}', fontsize=16, pad=25)
        plt.suptitle(f'Filtering threshold ($\\tau$): {threshold*100}% | Taxa shown: {len(keep_taxa)} + Others', fontsize=12, color='black', y=0.88)
        
        plt.ylabel('Relative Abundance', fontsize=12, fontweight='bold')
        plt.xlabel(f'Category: {category_col}', fontsize=12, fontweight='bold')
        plt.xticks(rotation=0, ha='center', fontsize=11) # Categories usually require less rotation
        
        plt.legend(title=rank_level.capitalize(), bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9, borderaxespad=0.)
        plt.tight_layout()
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Microbiome Relative Abundance Barplot Generator")
    parser.add_argument('-d', '--data', type=str, required=True, help='Path to taxa Excel file.')
    parser.add_argument('-r', '--rank', type=str, default='genus', help='Taxonomic rank (e.g., species, genus, family, order).')
    parser.add_argument('-t', '--threshold', type=float, default=0.01, help='Abundance threshold (default: 0.01).')
    parser.add_argument('-o', '--output', type=str, help='Custom output PDF filename.')
    
    # Grouping arguments
    parser.add_argument('-m', '--metadata', type=str, help='Path to metadata CSV file for grouping.')
    parser.add_argument('-c', '--category', type=str, help='Metadata column to group by (e.g., SEX, TYPE).')
    parser.add_argument('-id', '--sample_id', type=str, default='SampleID', help='Metadata column name containing sample IDs (default: SampleID).')

    args = parser.parse_args()
    if args.metadata and args.category:
        out_file = args.output if args.output else f"Microbiome_{args.rank}_{args.category}_{args.threshold}.pdf"
        print(f"Generating grouped barplot: {out_file}")
        generate_grouped_microbiome_pdf(
            data_path=args.data,
            metadata_path=args.metadata,
            category_col=args.category,
            sample_id_col=args.sample_id,
            rank_level=args.rank,
            threshold=args.threshold,
            output_file=out_file
        )
    else:
        out_file = args.output if args.output else f"Microbiome_{args.rank}_individual_{args.threshold}.pdf"
        print(f"Generating individual sample barplot: {out_file}")
        generate_microbiome_pdf(
            data_path=args.data,
            rank_level=args.rank,
            threshold=args.threshold,
            output_file=out_file
        )
    print("Execution complete.")

