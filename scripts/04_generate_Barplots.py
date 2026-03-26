import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

def export_topological_projection(fig, output_base, fmt, pad_inches=1.0):
    """
    Projects the continuous mathematical representation (matplotlib Figure) into
    a defined discrete or continuous state space (PNG, TIFF, or PDF).
    Maintains 300 DPI high-frequency spatial resolution for clinical rasters.
    """
    if fmt.lower() == 'pdf':
        with PdfPages(f"{output_base}.pdf") as pdf:
            pdf.savefig(fig, bbox_inches='tight', pad_inches=pad_inches, facecolor=fig.get_facecolor())
    elif fmt.lower() in ['png', 'tiff']:
        ext = 'tif' if fmt.lower() == 'tiff' else 'png'
        fig.savefig(f"{output_base}.{ext}", format=fmt.lower(), dpi=300, bbox_inches='tight', pad_inches=pad_inches, facecolor=fig.get_facecolor())
    else:
        raise ValueError(f"Unsupported topological projection format: {fmt}")

def generate_individual_microbiome_plot(data_path, rank_level='family', threshold=0.01, organism_name='Microbiome', output_base='microbiome_report', fmt='pdf'):
    """
    Generates a stacked bar chart of the relative abundance.
    Margins and 'Others' category are removed. Scaled precisely to 100%.
    Includes global padding and enlarged, square legend handles.
    """
    # 1. Load and segment data by Taxonomic Level
    df = pd.read_excel(data_path, sheet_name=0) 
 
    df_rank = df[df['Rank'].str.lower() == rank_level.lower()].copy()
    
    if df_rank.empty:
        raise ValueError(f"No data found for the taxonomic level: {rank_level}")
        
    df_rank.set_index('Scientific Name', inplace=True)
    
    sample_cols = [col for col in df_rank.columns if col not in ['Rank', 'TaxID', 'original_header']]
    df_counts = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # 2. Transformation to Relative Abundance
    rel_abund = df_counts.div(df_counts.sum(axis=0), axis=1)
    
    # 3. Mean calculation and Threshold Filtering
    mean_rel_abund = rel_abund.mean(axis=1)
    keep_taxa = mean_rel_abund[mean_rel_abund >= threshold].index
    
    df_kept = rel_abund.loc[keep_taxa]
    
    # --- RE-NORMALIZATION STEP ---
    df_kept_normalized = df_kept.div(df_kept.sum(axis=0), axis=1) * 100
    df_plot = df_kept_normalized.T
    
    colors = [
        "#b988d5", "#cbd588", "#88a05b", "#ffe156", "#6b8ba4", "#fda47d", 
        "#8e5634", "#d44645", "#5d9ca3", "#63b7af", "#dcd3ff", "#ff94cc", 
        "#ffa45b", "#806d40", "#2a363b", "#99b898", "#feceab", "#ff847c", 
        "#e84a5f", "#2a363b", "#56a5cc", "#c3e88d", "#ffcc5c", "#b09f59", 
        "#ff5e57", "#674d3c", "#4c4f69", "#8372a8", "#ff7c43", "#6a8a82"
    ]
    
    # 4. Figure Creation and Vector Rendering
    fig, ax = plt.subplots(figsize=(16, 12))
    
    plot_colors = colors[:len(keep_taxa)]
    
    df_plot.plot(kind='bar', stacked=True, ax=ax, width=0.85, color=plot_colors)
    ax.set_ylim(0, 100)
    
    for spine in ax.spines.values():
        spine.set_visible(False)
        
    ax.tick_params(axis='both', which='both', length=0)
    
    plt.ylabel('{} Relative Abundance >{:g}%'.format(organism_name, threshold * 100), fontsize=12, fontweight='bold')
    plt.xlabel('Sample Identifier', fontsize=12, fontweight='bold')
    
    plt.xticks(rotation=45, ha='right', fontsize=9)
    
    # Updated Legend
    plt.legend(title=rank_level.capitalize(), 
               bbox_to_anchor=(1.02, 1), 
               loc='upper left', 
               fontsize=11, 
               title_fontproperties={'weight': 'bold', 'size': 11},
               borderaxespad=0.,
               frameon=False,
               alignment='left',
               handlelength=1.5,
               handleheight=1.5,
               handletextpad=0.8)
    
    plt.tight_layout(pad=3.0)
    export_topological_projection(fig, output_base, fmt, pad_inches=1.0)
    plt.close(fig)


def generate_grouped_microbiome_plots(data_path, metadata_path, category_col, sample_id_col='SampleID', rank_level='family', threshold=0.01, organism_name='Microbiome', output_base='output', fmt='pdf'):
    """
    Groups absolute counts by a metadata category before calculating relative abundance.
    """
    df = pd.read_excel(data_path, sheet_name=0)
    
    if metadata_path.lower().endswith('.csv'):
        meta_df = pd.read_csv(metadata_path)
    elif metadata_path.lower().endswith(('.xls', '.xlsx')):
        meta_df = pd.read_excel(metadata_path)
    else:
        raise ValueError("Metadata file must be a .csv or .xlsx")
    
    df_rank = df[df['Rank'].str.lower() == rank_level.lower()].copy()
    if df_rank.empty:
        raise ValueError(f"No data found for the taxonomic level: {rank_level}")
        
    sample_cols = [col for col in df_rank.columns if col not in ['Rank', 'TaxID', 'original_header', 'Name']]
    df_rank[sample_cols] = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    df_counts = df_rank.groupby('Name')[sample_cols].sum()
    meta_df[sample_id_col] = meta_df[sample_id_col].astype(str).str.strip()
    category_map = dict(zip(meta_df[sample_id_col], meta_df[category_col]))
    
    df_counts.columns = [str(col).rsplit('_', 1)[0].strip() if '_' in str(col) else str(col).strip() for col in df_counts.columns]
    cols_to_keep = [col for col in df_counts.columns if col in category_map]
    
    if not cols_to_keep:
        raise ValueError("CRITICAL ERROR: None of the Sample IDs in your data matched the metadata.")
    
    df_counts = df_counts[cols_to_keep]
    df_counts_mapped = df_counts.rename(columns=category_map)
    df_grouped_counts = df_counts_mapped.groupby(df_counts_mapped.columns, axis=1).sum()
    
    rel_abund = df_grouped_counts.div(df_grouped_counts.sum(axis=0), axis=1)
    
    mean_rel_abund = rel_abund.mean(axis=1)
    keep_taxa = mean_rel_abund[mean_rel_abund >= threshold].index
    
    df_kept = rel_abund.loc[keep_taxa]
    
    df_kept_normalized = df_kept.div(df_kept.sum(axis=0), axis=1) * 100
    df_plot = df_kept_normalized.T
    
    colors = [ 
        "#b988d5", "#cbd588", "#88a05b", "#ffe156", "#6b8ba4",
        "#fda47d", "#8e5634", "#d44645", "#5d9ca3", "#63b7af",
        "#dcd3ff", "#ff94cc", "#ffa45b", "#806d40", "#2a363b",
        "#99b898", "#feceab", "#ff847c", "#e84a5f", "#2a363b",
        "#56a5cc", "#c3e88d", "#ffcc5c", "#b09f59", "#ff5e57",
        "#674d3c", "#4c4f69", "#8372a8", "#ff7c43", "#6a8a82",
        "#bfc0c0", "#de6e47", "#87bdd8", "#daedbd", "#ffb88c",
        "#fb8b24", "#4a7c59", "#ef476f", "#118ab2", "#06d6a0",
        "#f4a261", "#2a9d8f", "#264653", "#e63946", "#7d8491" 
    ]

    fig, ax = plt.subplots(figsize=(16, 12))
    plot_colors = colors[:len(keep_taxa)]
    
    df_plot.plot(kind='bar', stacked=True, ax=ax, width=0.5, color=plot_colors)
    
    ax.set_ylim(0, 100)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['left'].set_color('black')
    ax.spines['left'].set_linewidth(1.7)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_color('black')
    ax.spines['bottom'].set_linewidth(1.7)
    ax.spines['bottom'].set_position(('outward', 25))

    ax.tick_params(axis='y', which='both', length=8, color='black')
    ax.tick_params(axis='x', which='both', length=8, color='black')
    
    plt.setp(ax.get_yticklabels(), fontweight='bold', fontsize=11)
    plt.ylabel('{} relative abundance (>{:g}%)'.format(organism_name, threshold * 100), fontsize=17, fontweight='bold')
    plt.xlabel(f'{category_col}', fontsize=17, fontweight='bold')
    plt.xticks(rotation=0, ha='center', fontsize=14)
    
    plt.legend(title=rank_level.capitalize(), 
               bbox_to_anchor=(1., .80), 
               loc='upper left', 
               prop={'style': 'italic', 'size': 12},
               title_fontproperties={'weight': 'bold', 'size': 15}, 
               borderaxespad=0.,
               frameon=False,
               alignment='left',
               handlelength=1.5,
               handleheight=1.5,
               handletextpad=0.8)
               
    plt.tight_layout(pad=3.0)
    export_topological_projection(fig, output_base, fmt, pad_inches=1.0)
    plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Relative Abundance Barplot Generator")
    parser.add_argument('-d', '--data', type=str, required=True, help='Path to taxa Excel file.')
    parser.add_argument('-r', '--rank', type=str, default='genus', help='Taxonomic rank.')
    parser.add_argument('-t', '--threshold', type=float, default=0.01, help='Abundance threshold.')
    parser.add_argument('-o', '--output', type=str, help='Custom output filename base.')
    parser.add_argument('-org', '--organism', type=str, default='Microbiome', help='Prefix for the Y-axis label.')
    parser.add_argument('-m', '--metadata', type=str, help='Path to metadata CSV file for grouping.')
    parser.add_argument('-c', '--category', type=str, help='Metadata column to group by.')
    parser.add_argument('-id', '--sample_id', type=str, default='SampleID', help='Metadata column name for sample IDs.')
    
    parser.add_argument('-fmt', '--format', type=str, choices=['pdf', 'png', 'tiff'], default='pdf', 
                        help='Topological projection format: pdf (continuous vector), png/tiff (discrete raster at 300 DPI).')

    args = parser.parse_args()
    
    safe_org_name = args.organism.replace(" ", "_")
    output_base = args.output if args.output else f"{safe_org_name}_{args.rank}_{args.threshold}"
    
    if args.metadata and args.category:
        generate_grouped_microbiome_plots(
            data_path=args.data,
            metadata_path=args.metadata,
            category_col=args.category,
            sample_id_col=args.sample_id,
            rank_level=args.rank,
            threshold=args.threshold,
            organism_name=args.organism, 
            output_base=output_base,
            fmt=args.format.lower()
        )
    else:
        generate_individual_microbiome_plot(
            data_path=args.data,
            rank_level=args.rank,
            threshold=args.threshold,
            organism_name=args.organism, 
            output_base=output_base,
            fmt=args.format.lower()
        )
