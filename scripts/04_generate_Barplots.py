import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

# ---------------------------------------------------------------------------
# Unified color palette — identical across all pipeline scripts (05, 06, 07)
# ---------------------------------------------------------------------------
COLORS = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#b988d5", "#cbd588", "#88a05b", "#ffe156", "#6b8ba4", "#fda47d",
    "#8e5634", "#d44645", "#5d9ca3", "#63b7af", "#dcd3ff", "#ff94cc",
    "#ffa45b", "#806d40", "#2a363b", "#99b898", "#feceab", "#ff847c",
    "#e84a5f", "#2a363b", "#56a5cc", "#c3e88d", "#ffcc5c", "#b09f59",
    "#ff5e57", "#674d3c", "#4c4f69", "#8372a8", "#ff7c43", "#6a8a82"
]


def export_topological_projection(fig, output_base, fmt, pad_inches=1.0):
    """
    Projects the continuous mathematical representation (matplotlib Figure) into
    a defined discrete or continuous state space (PNG, TIFF, or PDF).
    Maintains 300 DPI high-frequency spatial resolution for clinical rasters.
    """
    if fmt.lower() == 'pdf':
        with PdfPages(f"{output_base}.pdf") as pdf:
            pdf.savefig(fig, bbox_inches='tight', pad_inches=pad_inches,
                        facecolor=fig.get_facecolor())
    elif fmt.lower() in ['png', 'tiff']:
        ext = 'tif' if fmt.lower() == 'tiff' else 'png'
        fig.savefig(f"{output_base}.{ext}", format=fmt.lower(), dpi=300,
                    bbox_inches='tight', pad_inches=pad_inches,
                    facecolor=fig.get_facecolor())
    else:
        raise ValueError(f"Unsupported topological projection format: {fmt}")


def export_summary_table(df, output_base, sheet_name='Abundance_Summary'):
    """
    Serializes the computed abundance summary DataFrame into a formatted
    Excel workbook with auto-adjusted column widths for readability.
    """
    output_path = f"{output_base}_table.xlsx"
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
        ws = writer.sheets[sheet_name]
        for col in ws.columns:
            max_len = max(
                (len(str(cell.value)) if cell.value is not None else 0)
                for cell in col
            )
            ws.column_dimensions[col[0].column_letter].width = min(max_len + 4, 50)
    print(f"    [*] Summary table saved: {output_path}")


def generate_individual_microbiome_plot(data_path, rank_level='family', threshold=0.01,
                                         organism_name='Microbiome',
                                         output_base='microbiome_report',
                                         fmt='pdf', sample_order=None, no_table=False):
    """
    Generates a stacked bar chart of the relative abundance per individual sample.
    Margins and 'Others' category are removed. Scaled precisely to 100%.
    Includes global padding and enlarged, square legend handles.
    Optionally exports an abundance summary table (.xlsx).
    """
    # 1. Load and segment data by Taxonomic Level
    df = pd.read_excel(data_path, sheet_name=0)
    df_rank = df[df['Rank'].str.lower() == rank_level.lower()].copy()

    if df_rank.empty:
        raise ValueError(f"No data found for the taxonomic level: {rank_level}")

    df_rank.set_index('Scientific Name', inplace=True)

    sample_cols = [col for col in df_rank.columns
                   if col not in ['Rank', 'TaxID', 'original_header']]
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

    # --- CUSTOM INDEX REORDERING ---
    if sample_order:
        valid_order = [item for item in sample_order if item in df_plot.index]
        missing     = [item for item in df_plot.index if item not in valid_order]
        df_plot     = df_plot.loc[valid_order + missing]

    # 4. Figure Creation and Vector Rendering
    fig, ax = plt.subplots(figsize=(16, 12))

    plot_colors = COLORS[:len(keep_taxa)]

    df_plot.plot(kind='bar', stacked=True, ax=ax, width=0.85, color=plot_colors)
    ax.set_ylim(0, 100)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.tick_params(axis='both', which='both', length=0)

    plt.ylabel('{} Relative Abundance >{:g}%'.format(organism_name, threshold * 100),
               fontsize=12, fontweight='bold')
    plt.xlabel('Sample Identifier', fontsize=12, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=9)

    # Determine if legend items should be italicized based on biological taxonomy rank rules
    font_style = 'italic' if rank_level.lower() in ['genus', 'species'] else 'normal'

    plt.legend(
        title=rank_level.capitalize(),
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        prop={'style': font_style, 'size': 11},
        title_fontproperties={'weight': 'bold', 'size': 11},
        borderaxespad=0.,
        frameon=False,
        alignment='left',
        handlelength=1.5,
        handleheight=1.5,
        handletextpad=0.8
    )

    plt.tight_layout(pad=3.0)
    export_topological_projection(fig, output_base, fmt, pad_inches=1.0)
    plt.close(fig)

    # 5. Abundance Summary Table — Taxon | Mean_Rel_Abundance_pct | per-sample columns
    if not no_table:
        summary_df = df_plot.T.copy()
        summary_df.insert(0, 'Taxon', summary_df.index)
        summary_df.insert(1, 'Mean_Rel_Abundance_pct',
                          df_plot.mean(axis=0).values.round(4))
        summary_df = (
            summary_df
            .round(4)
            .sort_values('Mean_Rel_Abundance_pct', ascending=False)
            .reset_index(drop=True)
        )
        export_summary_table(summary_df, output_base)


def generate_grouped_microbiome_plots(data_path, metadata_path, category_col,
                                       sample_id_col='SampleID', rank_level='family',
                                       threshold=0.01, organism_name='Microbiome',
                                       output_base='output', fmt='pdf',
                                       category_order=None, no_table=False):
    """
    Groups absolute counts by a metadata category before calculating relative abundance.
    Includes an optional custom categorical order mapping.
    Optionally exports an abundance summary table (.xlsx).
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

    sample_cols = [col for col in df_rank.columns
                   if col not in ['Rank', 'TaxID', 'original_header', 'Name']]
    df_rank[sample_cols] = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

    df_counts    = df_rank.groupby('Name')[sample_cols].sum()
    meta_df[sample_id_col] = meta_df[sample_id_col].astype(str).str.strip()
    category_map = dict(zip(meta_df[sample_id_col], meta_df[category_col]))

    df_counts.columns = [
        str(col).rsplit('_', 1)[0].strip() if '_' in str(col) else str(col).strip()
        for col in df_counts.columns
    ]
    cols_to_keep = [col for col in df_counts.columns if col in category_map]

    if not cols_to_keep:
        raise ValueError("CRITICAL ERROR: None of the Sample IDs in your data matched the metadata.")

    df_counts        = df_counts[cols_to_keep]
    df_counts_mapped = df_counts.rename(columns=category_map)

    # --- Modern pandas equivalent of .groupby(columns, axis=1).sum() ---
    df_grouped_counts = df_counts_mapped.T.groupby(df_counts_mapped.T.index).sum().T

    rel_abund      = df_grouped_counts.div(df_grouped_counts.sum(axis=0), axis=1)
    mean_rel_abund = rel_abund.mean(axis=1)
    keep_taxa      = mean_rel_abund[mean_rel_abund >= threshold].index

    df_kept            = rel_abund.loc[keep_taxa]
    df_kept_normalized = df_kept.div(df_kept.sum(axis=0), axis=1) * 100
    df_plot            = df_kept_normalized.T

    # --- CUSTOM INDEX REORDERING ---
    if category_order:
        valid_order = [cat for cat in category_order if cat in df_plot.index]
        missing     = [cat for cat in df_plot.index if cat not in valid_order]
        df_plot     = df_plot.loc[valid_order + missing]

    fig, ax     = plt.subplots(figsize=(16, 12))
    plot_colors = COLORS[:len(keep_taxa)]

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
    plt.ylabel('{} relative abundance (>{:g}%)'.format(organism_name, threshold * 100),
               fontsize=17, fontweight='bold')
    plt.xlabel(f'{category_col}', fontsize=17, fontweight='bold')
    plt.xticks(rotation=0, ha='center', fontsize=14)

    # Determine if legend items should be italicized based on biological taxonomy rank rules
    font_style = 'italic' if rank_level.lower() in ['genus', 'species'] else 'normal'

    plt.legend(
        title=rank_level.capitalize(),
        bbox_to_anchor=(1., .80),
        loc='upper left',
        prop={'style': font_style, 'size': 12},
        title_fontproperties={'weight': 'bold', 'size': 15},
        borderaxespad=0.,
        frameon=False,
        alignment='left',
        handlelength=1.5,
        handleheight=1.5,
        handletextpad=0.8
    )

    plt.tight_layout(pad=3.0)
    export_topological_projection(fig, output_base, fmt, pad_inches=1.0)
    plt.close(fig)

    # Abundance Summary Table — Taxon | Mean_Rel_Abundance_pct | per-group columns
    if not no_table:
        summary_df = df_plot.T.copy()
        summary_df.insert(0, 'Taxon', summary_df.index)
        summary_df.insert(1, 'Mean_Rel_Abundance_pct',
                          df_plot.mean(axis=0).values.round(4))
        summary_df = (
            summary_df
            .round(4)
            .sort_values('Mean_Rel_Abundance_pct', ascending=False)
            .reset_index(drop=True)
        )
        export_summary_table(summary_df, output_base)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Relative Abundance Barplot Generator")
    parser.add_argument('-d',    '--data',      type=str, required=True,
                        help='Path to taxa Excel file.')
    parser.add_argument('-r',    '--rank',      type=str, default='genus',
                        help='Taxonomic rank.')
    parser.add_argument('-t',    '--threshold', type=float, default=0.01,
                        help='Abundance threshold (e.g. 0.01 = 1%%).')
    parser.add_argument('-o',    '--output',    type=str,
                        help='Custom output filename base (no extension).')
    parser.add_argument('-org',  '--organism',  type=str, default='Microbiome',
                        help='Prefix for the Y-axis label.')
    parser.add_argument('-m',    '--metadata',  type=str,
                        help='Path to metadata CSV or Excel file for grouping.')
    parser.add_argument('-c',    '--category',  type=str,
                        help='Metadata column to group by.')
    parser.add_argument('-id',   '--sample_id', type=str, default='SampleID',
                        help='Metadata column name for sample IDs.')
    parser.add_argument('-fmt',  '--format',    type=str,
                        choices=['pdf', 'png', 'tiff'], default='pdf',
                        help='Output format: pdf (vector), png/tiff (raster 300 DPI).')
    parser.add_argument('-ord',  '--order',     type=str, nargs='+',
                        help='Custom order for x-axis categories (space-separated).')
    parser.add_argument('--no_table', action='store_true',
                        help='Skip exporting the abundance summary table (.xlsx).')

    args = parser.parse_args()

    OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'results', 'Barplots')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    safe_org_name = args.organism.replace(" ", "_")
    base_name     = args.output if args.output else f"{safe_org_name}_{args.rank}_{args.threshold}"
    output_base   = os.path.join(OUTPUT_DIR, base_name)

    if args.metadata and args.category:
        generate_grouped_microbiome_plots(
            data_path     = args.data,
            metadata_path = args.metadata,
            category_col  = args.category,
            sample_id_col = args.sample_id,
            rank_level    = args.rank,
            threshold     = args.threshold,
            organism_name = args.organism,
            output_base   = output_base,
            fmt           = args.format.lower(),
            category_order= args.order,
            no_table      = args.no_table,
        )
    else:
        generate_individual_microbiome_plot(
            data_path     = args.data,
            rank_level    = args.rank,
            threshold     = args.threshold,
            organism_name = args.organism,
            output_base   = output_base,
            fmt           = args.format.lower(),
            sample_order  = args.order,
            no_table      = args.no_table,
        )
