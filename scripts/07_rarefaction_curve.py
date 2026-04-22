import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import math

# ---------------------------------------------------------------------------
# Unified color palette — identical across all pipeline scripts
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

# ---------------------------------------------------------------------------
# Global typography
# ---------------------------------------------------------------------------
matplotlib.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.weight":      "bold",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "axes.labelsize":   11,
    "axes.titlesize":   12,
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
})

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def apply_theme_bw(ax):
    """White panel background, 4 black spines, outward ticks."""
    ax.set_facecolor("white")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")
    ax.tick_params(axis="both", which="both", direction="out",
                   colors="black", width=0.8, length=4)
    ax.tick_params(axis="x", labelcolor="black", labelsize=9)
    ax.tick_params(axis="y", labelcolor="black", labelsize=9)

def export_figure(fig, output_base, fmt, pad_inches=0.35):
    """Exports matplotlib figure to PDF, PNG, or TIFF."""
    if fmt == "pdf":
        with PdfPages(f"{output_base}.pdf") as pdf:
            pdf.savefig(fig, bbox_inches="tight", pad_inches=pad_inches,
                        facecolor=fig.get_facecolor())
    elif fmt in ("png", "tiff"):
        ext = "tif" if fmt == "tiff" else "png"
        fig.savefig(f"{output_base}.{ext}", format=fmt, dpi=300,
                    bbox_inches="tight", pad_inches=pad_inches,
                    facecolor=fig.get_facecolor())
    else:
        raise ValueError(f"Unsupported format: {fmt}")

def export_summary_table(df, output_base, sheet_name='Summary'):
    """Serializes DataFrame to auto-formatted Excel."""
    output_path = f"{output_base}_table.xlsx"
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
        ws = writer.sheets[sheet_name]
        for col in ws.columns:
            max_len = max((len(str(cell.value)) if cell.value is not None else 0) for cell in col)
            ws.column_dimensions[col[0].column_letter].width = min(max_len + 4, 50)
    print(f"    [*] Summary table saved: {output_path}")

# ---------------------------------------------------------------------------
# Rarefaction Math & Simulation
# ---------------------------------------------------------------------------
def simulate_rarefaction(counts_series, depth_steps, n_iter):
    """
    Subsamples reads iteratively without replacement using np.random.choice.
    Returns means and standard deviations of observed species richness.
    """
    counts = counts_series.values.astype(int)
    counts = counts[counts > 0]
    total_reads = counts.sum()

    if total_reads == 0:
        return np.zeros(len(depth_steps)), np.zeros(len(depth_steps))

    # Construct the population array of reads: e.g. [taxon0, taxon0, taxon1...]
    taxa_indices = np.arange(len(counts))
    reads_population = np.repeat(taxa_indices, counts)

    means = []
    stds = []

    for depth in depth_steps:
        depth = int(depth)
        if depth == 0:
            means.append(0.0)
            stds.append(0.0)
            continue
        
        if depth >= total_reads:
            means.append(float(len(counts)))
            stds.append(0.0)
            continue

        iter_richness = []
        for _ in range(n_iter):
            sampled = np.random.choice(reads_population, size=depth, replace=False)
            iter_richness.append(len(np.unique(sampled)))
        
        means.append(np.mean(iter_richness))
        stds.append(np.std(iter_richness, ddof=1) if n_iter > 1 else 0.0)

    return np.array(means), np.array(stds)

# ---------------------------------------------------------------------------
# Main Pipeline
# ---------------------------------------------------------------------------
def generate_rarefaction_curves(data_path, metadata_path, category_col, sample_id_col,
                                rank_level, target_depth, num_steps, n_iter,
                                organism_name, output_base, fmt, no_table):
    
    print(f"[*] Generating empirical rarefaction curves → {output_base}")

    # 1. Load Data
    df = pd.read_excel(data_path, sheet_name=0)
    df_rank = df[df["Rank"].str.lower() == rank_level.lower()].copy()
    if df_rank.empty:
        raise ValueError(f"No data found for taxonomic level: {rank_level}")

    meta_cols = {"Rank", "TaxID", "original_header", "Name", "Scientific Name"}
    sample_cols = [c for c in df_rank.columns if c not in meta_cols]
    df_counts = df_rank[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)

    # 2. Extract total reads per sample
    reads_totals = df_counts.sum(axis=0)
    min_total_reads = int(reads_totals.min())
    
    if target_depth is None:
        target_depth = min_total_reads
        print(f"    [*] Auto-detected rarefaction depth (min reads): {target_depth}")
    else:
        target_depth = int(target_depth)
        print(f"    [*] User-supplied rarefaction depth: {target_depth}")

    # 3. Load Metadata (if provided)
    meta_map = {}
    use_groups = False
    
    if metadata_path and category_col:
        meta_df = (pd.read_csv(metadata_path) if metadata_path.lower().endswith(".csv") 
                   else pd.read_excel(metadata_path))
        meta_df[sample_id_col] = meta_df[sample_id_col].astype(str).str.strip().str.lower()

        if pd.api.types.is_numeric_dtype(meta_df[category_col]):
            meta_df[category_col] = pd.qcut(meta_df[category_col].dropna(), q=4).astype(str)

        meta_map = dict(zip(meta_df[sample_id_col], meta_df[category_col].astype(str)))
        use_groups = True

    # 4. Color Assignment
    if use_groups:
        groups = sorted(list(set(meta_map.values())))
        palette_map = {g: COLORS[i % len(COLORS)] for i, g in enumerate(groups)}
    else:
        palette_map = {s: COLORS[i % len(COLORS)] for i, s in enumerate(sample_cols)}

    # 5. Initialize Figure
    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
    apply_theme_bw(ax)

    summary_records = []
    
    # 6. Process each sample
    print(f"    [*] Computing {n_iter} iterations per depth for {len(sample_cols)} sample(s)...")
    for sample in sample_cols:
        clean_sample = str(sample).split("_")[0].strip().lower()
        group = meta_map.get(clean_sample, "Unknown") if use_groups else "N/A"
        color = palette_map.get(group) if use_groups else palette_map[sample]
        
        N = int(reads_totals[sample])
        S_obs = int((df_counts[sample] > 0).sum())

        # Define evaluation depths: uniform steps up to N, explicitly including target_depth
        base_steps = np.linspace(0, N, num=num_steps, dtype=int)
        eval_depths = np.unique(np.sort(np.append(base_steps, [target_depth])))
        eval_depths = eval_depths[eval_depths <= N] # Don't simulate past total reads
        
        means, stds = simulate_rarefaction(df_counts[sample], eval_depths, n_iter)
        
        # Plot Curve and Confidence Band
        ax.plot(eval_depths, means, color=color, linewidth=1.5, alpha=0.85)
        ax.fill_between(eval_depths, means - stds, means + stds, color=color, alpha=0.15, edgecolor='none')
        
        # Extract specific richness at the target rarefaction depth
        rich_at_depth = np.nan
        if target_depth <= N:
            idx = np.where(eval_depths == target_depth)[0][0]
            rich_at_depth = means[idx]

        summary_records.append({
            'Sample': sample,
            'Group': group,
            'Reads_Total': N,
            'Species_Observed': S_obs,
            f'Rarefied_Richness_at_{target_depth}': rich_at_depth
        })

    # 7. Formatting and vertical line
    ax.axvline(x=target_depth, color="#d32f2f", linestyle="--", linewidth=1.5, zorder=5, 
               label=f"Depth: {target_depth}")

    ax.set_xlabel("Number of Sequences Sampled (Depth)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Species Richness (Observed Taxa)", fontsize=11, fontweight="bold")
    ax.set_title(f"{organism_name} — Rarefaction Curves ({rank_level.capitalize()})", 
                 fontsize=13, fontweight="bold", pad=12)

    # 8. Legend
    if use_groups:
        legend_patches = [mpatches.Patch(color=palette_map[g], label=g) for g in groups]
        legend_patches.append(plt.Line2D([0], [0], color="#d32f2f", linestyle="--", label=f"Depth: {target_depth}"))
        ax.legend(handles=legend_patches, title=category_col, frameon=True, edgecolor="black")

    fig.tight_layout()
    
    # 9. Export Figure
    export_figure(fig, output_base, fmt)
    plt.close(fig)
    print(f"    [✓] Rarefaction plot saved: {output_base}.{fmt}")

    # 10. Summary Table
    if not no_table:
        table_df = pd.DataFrame(summary_records)
        export_summary_table(table_df, output_base, sheet_name='Rarefaction_Summary')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rarefaction Curve Generator (Stochastic Subsampling)")
    
    parser.add_argument("-d",    "--data",      required=True, help="Path to taxa Excel file.")
    parser.add_argument("-r",    "--rank",      default="genus", help="Taxonomic rank level.")
    parser.add_argument("-o",    "--output",    default=None, help="Output filename base.")
    parser.add_argument("-org",  "--organism",  default="Microbiome", help="Organism prefix for title.")
    
    parser.add_argument("-m",    "--metadata",  default=None, help="Path to metadata CSV/Excel (optional).")
    parser.add_argument("-c",    "--category",  default=None, help="Metadata column to group by (optional).")
    parser.add_argument("-id",   "--sample_id", default="SampleID", help="Metadata column containing sample IDs.")
    
    parser.add_argument("--depth", type=int,    default=None, help="Target rarefaction depth (default: min total reads).")
    parser.add_argument("--steps", type=int,    default=10, help="Number of depth intervals to compute (default: 10).")
    parser.add_argument("--iter",  type=int,    default=10, help="Monte Carlo resampling iterations per depth (default: 10).")
    
    parser.add_argument("-fmt",  "--format",    choices=["pdf", "png", "tiff"], default="pdf", help="Output format.")
    parser.add_argument("--no_table", action="store_true", help="Skip exporting summary table.")

    args = parser.parse_args()

    OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'Rarefaction')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    safe_name = args.organism.replace(" ", "_")
    base_name = args.output if args.output else f"{safe_name}_{args.rank}_Rarefaction"
    output_base = os.path.join(OUTPUT_DIR, base_name)

    if (args.metadata and not args.category) or (args.category and not args.metadata):
        print("[!] Warning: Both --metadata and --category must be provided for group coloring. Falling back to sample colors.")

    generate_rarefaction_curves(
        data_path     = args.data,
        metadata_path = args.metadata,
        category_col  = args.category,
        sample_id_col = args.sample_id,
        rank_level    = args.rank,
        target_depth  = args.depth,
        num_steps     = args.steps,
        n_iter        = args.iter,
        organism_name = args.organism,
        output_base   = output_base,
        fmt           = args.format.lower(),
        no_table      = args.no_table,
    )
