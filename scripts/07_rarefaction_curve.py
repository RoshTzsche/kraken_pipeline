import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

# ---------------------------------------------------------------------------
# Unified color palette — distinct colors for each sample
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
# Global typography settings for the plot
# ---------------------------------------------------------------------------
matplotlib.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.weight":      "bold",
    "axes.labelweight": "bold",
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "axes.labelsize":   11,
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
})

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def apply_theme_minimal(ax):
    """Applies a highly minimal theme to the matplotlib axes."""
    ax.set_facecolor("white")
    ax.grid(False)
    
    # Remove top and right spines for a minimal look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Keep bottom and left spines
    ax.spines['bottom'].set_linewidth(0.8)
    ax.spines['bottom'].set_color("black")
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['left'].set_color("black")
    
    ax.tick_params(axis="both", which="both", direction="out",
                   colors="black", width=0.8, length=4)

def export_figure(fig, output_base, fmt, pad_inches=0.1):
    """Exports the matplotlib figure to the specified format (PDF, PNG, TIFF)."""
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
    """Serializes the summary DataFrame to an auto-formatted Excel file."""
    output_path = f"{output_base}_table.xlsx"
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
        ws = writer.sheets[sheet_name]
        # Auto-adjust column widths
        for col in ws.columns:
            max_len = max((len(str(cell.value)) if cell.value is not None else 0) for cell in col)
            ws.column_dimensions[col[0].column_letter].width = min(max_len + 4, 50)
    print(f"    [*] Summary table saved: {output_path}")

# ---------------------------------------------------------------------------
# Rarefaction Math & Simulation
# ---------------------------------------------------------------------------
def simulate_rarefaction(counts_series, depth_steps, n_iter):
    """
    Subsamples reads iteratively without replacement.
    Returns means and standard deviations of observed species richness.
    """
    counts = counts_series.values.astype(int)
    counts = counts[counts > 0]
    total_reads = counts.sum()

    if total_reads == 0:
        return np.zeros(len(depth_steps)), np.zeros(len(depth_steps))

    # Construct the population array of reads
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
        
        # If requested depth exceeds total reads, it caps at total observed richness
        if depth >= total_reads:
            means.append(float(len(counts)))
            stds.append(0.0)
            continue

        # Monte Carlo resampling
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
def generate_rarefaction_curves(data_path, rank_level, target_depth, num_steps, n_iter,
                                organism_name, output_base, fmt, no_table):
    
    print(f"[*] Generating empirical rarefaction curves (Global Mode) → {output_base}")

    # 1. Load Data
    df = pd.read_excel(data_path, sheet_name=0)
    df_rank = df[df["Rank"].str.lower() == rank_level.lower()].copy()
    if df_rank.empty:
        raise ValueError(f"No data found for taxonomic level: {rank_level}")

    # Separate metadata columns from sample columns
    meta_cols = {"Rank", "TaxID", "original_header", "Name", "Scientific Name"}
    sample_cols = [c for c in df_rank.columns if c not in meta_cols]
    df_counts = df_rank[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)

    # 2. Extract total reads per sample
    reads_totals = df_counts.sum(axis=0)
    min_total_reads = int(reads_totals.min())
    
    # Auto-detect target depth if not provided
    if target_depth is None:
        target_depth = min_total_reads
        print(f"    [*] Auto-detected rarefaction depth (min reads): {target_depth}")
    else:
        target_depth = int(target_depth)
        print(f"    [*] User-supplied rarefaction depth: {target_depth}")

    # 3. Color Assignment (Strictly one unique color per sample)
    palette_map = {sample: COLORS[i % len(COLORS)] for i, sample in enumerate(sample_cols)}

    # 4. Initialize Figure
    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
    apply_theme_minimal(ax)

    summary_records = []
    
    # 5. Process each sample individually
    print(f"    [*] Computing {n_iter} iterations per depth for {len(sample_cols)} sample(s)...")
    for sample in sample_cols:
        color = palette_map[sample]
        
        N = int(reads_totals[sample])
        S_obs = int((df_counts[sample] > 0).sum())

        # Define evaluation depths: 
        # High resolution (num_steps) to ensure the curve smoothly visually flattens (plateaus)
        base_steps = np.linspace(0, N, num=num_steps, dtype=int)
        
        # Include the target depth explicitly if it's within the range
        eval_depths = np.unique(np.sort(np.append(base_steps, [target_depth])))
        eval_depths = eval_depths[eval_depths <= N] 
        
        # Ensure we always plot the absolute maximum to show the flattening effect completely
        if N not in eval_depths:
            eval_depths = np.append(eval_depths, N)
            
        means, stds = simulate_rarefaction(df_counts[sample], eval_depths, n_iter)
        
        # Plot the main curve and the confidence band
        ax.plot(eval_depths, means, color=color, linewidth=1.5, alpha=0.85)
        ax.fill_between(eval_depths, means - stds, means + stds, color=color, alpha=0.15, edgecolor='none')
        
        # Extract specific richness at the target rarefaction depth
        rich_at_depth = np.nan
        if target_depth <= N:
            # Find the closest depth in our evaluation array to the target
            idx = np.abs(eval_depths - target_depth).argmin()
            rich_at_depth = means[idx]

        summary_records.append({
            'Sample': sample,
            'Reads_Total': N,
            'Species_Observed': S_obs,
            f'Rarefied_Richness_at_{target_depth}': rich_at_depth
        })

    # 6. Formatting minimal axes and vertical depth line
    # Minimal indicator for target depth
    ax.axvline(x=target_depth, color="gray", linestyle=":", linewidth=1.2, zorder=1)
    ax.text(target_depth + (ax.get_xlim()[1] * 0.01), ax.get_ylim()[0], f"{target_depth}", 
            color="gray", fontsize=9, va='bottom', ha='left')

    # Strictly set axes names as requested
    ax.set_xlabel("depth", fontsize=12, fontweight="bold")
    ax.set_ylabel("observed taxa", fontsize=12, fontweight="bold")
    
    # 7. Export Figure (No title, no legends)
    export_figure(fig, output_base, fmt)
    plt.close(fig)
    print(f"    [✓] Rarefaction plot saved: {output_base}.{fmt}")

    # 8. Summary Table
    if not no_table:
        table_df = pd.DataFrame(summary_records)
        export_summary_table(table_df, output_base, sheet_name='Rarefaction_Summary')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Global Rarefaction Curve Generator (Stochastic Subsampling)")
    
    parser.add_argument("-d",    "--data",      required=True, help="Path to taxa Excel file.")
    parser.add_argument("-r",    "--rank",      default="genus", help="Taxonomic rank level.")
    parser.add_argument("-o",    "--output",    default=None, help="Output filename base.")
    parser.add_argument("-org",  "--organism",  default="Microbiome", help="Organism prefix for title.")
    
    parser.add_argument("--depth", type=int,    default=None, help="Target rarefaction depth (default: min total reads).")
    
    # Increased default steps to 50 to guarantee smooth curves that visibly flatten out
    parser.add_argument("--steps", type=int,    default=50, help="Number of depth intervals to compute (default: 50).")
    parser.add_argument("--iter",  type=int,    default=10, help="Monte Carlo resampling iterations per depth (default: 10).")
    
    parser.add_argument("-fmt",  "--format",    choices=["pdf", "png", "tiff"], default="pdf", help="Output format.")
    parser.add_argument("--no_table", action="store_true", help="Skip exporting summary table.")

    args = parser.parse_args()

    # Create output directory
    OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'Rarefaction')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    safe_name = args.organism.replace(" ", "_")
    base_name = args.output if args.output else f"{safe_name}_{args.rank}_Rarefaction_Global"
    output_base = os.path.join(OUTPUT_DIR, base_name)

    generate_rarefaction_curves(
        data_path     = args.data,
        rank_level    = args.rank,
        target_depth  = args.depth,
        num_steps     = args.steps,
        n_iter        = args.iter,
        organism_name = args.organism,
        output_base   = output_base,
        fmt           = args.format.lower(),
        no_table      = args.no_table,
    )