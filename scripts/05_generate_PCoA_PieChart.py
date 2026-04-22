import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.distance import pdist, squareform
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import scipy.stats as stats

# ---------------------------------------------------------------------------
# Unified color palette — identical across all pipeline scripts (04, 06, 07)
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


def export_topological_projection(fig, output_base, fmt, pad_inches=0.5):
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


def export_summary_table(df, output_base, sheet_name='Summary', extra_sheets=None):
    """
    Serializes one or more computed DataFrames into a formatted Excel workbook.
    extra_sheets: dict of {sheet_name: DataFrame} for additional worksheets
                  (e.g., Bray-Curtis distance matrix as a second sheet).
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

        if extra_sheets:
            for sname, sdf in extra_sheets.items():
                sname_safe = sname[:31]   # Excel worksheet name limit
                sdf.to_excel(writer, index=False, sheet_name=sname_safe)
                ws2 = writer.sheets[sname_safe]
                for col in ws2.columns:
                    max_len = max(
                        (len(str(cell.value)) if cell.value is not None else 0)
                        for cell in col
                    )
                    ws2.column_dimensions[col[0].column_letter].width = min(max_len + 4, 50)
    print(f"    [*] Summary table saved: {output_path}")


def compute_anosim(dist_matrix, groups, n_perm=999, seed=42):
    """
    Analysis of Similarities (ANOSIM) — Clarke 1993.
    R = (r̄_B − r̄_W) / (N(N−1)/4)
      where r̄_B = mean rank of between-group distances,
            r̄_W = mean rank of within-group distances.
    Permutation-based p-value (n_perm label shuffles).
    Returns: (R_statistic, p_value)
    """
    np.random.seed(seed)
    n          = len(groups)
    groups_arr = np.array(groups)

    if len(np.unique(groups_arr)) < 2:
        print("    [!] ANOSIM: fewer than 2 unique groups — statistic cannot be computed.")
        return np.nan, np.nan

    # Rank all condensed upper-triangle distances, then rebuild square rank matrix
    dist_flat     = squareform(dist_matrix, checks=False)
    ranked_flat   = stats.rankdata(dist_flat)
    rank_matrix   = squareform(ranked_flat, checks=False)
    triu_idx      = np.triu_indices(n, k=1)

    def _r(grps):
        same_mask = (grps[:, None] == grps[None, :])[triu_idx]
        r_triu    = rank_matrix[triu_idx]
        r_W = np.mean(r_triu[same_mask])   if same_mask.any()  else 0.0
        r_B = np.mean(r_triu[~same_mask])  if (~same_mask).any() else 0.0
        return (r_B - r_W) / (n * (n - 1) / 4)

    obs_r  = _r(groups_arr)
    perm_r = np.array([_r(np.random.permutation(groups_arr)) for _ in range(n_perm)])
    p_val  = (np.sum(perm_r >= obs_r) + 1) / (n_perm + 1)

    print(f"    [*] ANOSIM     R = {obs_r:.4f}   p = {p_val:.4f}   (n_perm = {n_perm})")
    return float(obs_r), float(p_val)


def compute_permanova(dist_matrix, groups, n_perm=999, seed=42):
    """
    Permutational MANOVA (PERMANOVA / adonis) — Anderson 2001.
    Pseudo-F = (SS_A / df_A) / (SS_W / df_W)
      where SS is computed on squared Bray-Curtis distances,
            df_A = a − 1  (between-group df),
            df_W = N − a  (within-group  df).
    Permutation-based p-value (n_perm label shuffles).
    Returns: (F_statistic, R2, p_value)
    """
    np.random.seed(seed)
    n           = len(groups)
    groups_arr  = np.array(groups)
    a           = len(np.unique(groups_arr))

    if a < 2 or n <= a:
        print("    [!] PERMANOVA: insufficient groups or samples — statistic cannot be computed.")
        return np.nan, np.nan, np.nan

    D    = dist_matrix ** 2       # squared distance matrix
    df_a = a - 1
    df_w = n - a

    def _pseudo_f(grps):
        ss_t = np.sum(D) / (2 * n)
        ss_w = 0.0
        for g in np.unique(grps):
            idx = np.where(grps == g)[0]
            n_g = len(idx)
            if n_g > 1:
                ss_w += np.sum(D[np.ix_(idx, idx)]) / (2 * n_g)
        ss_a = ss_t - ss_w
        if df_w <= 0 or ss_w == 0:
            return np.nan, np.nan
        return (ss_a / df_a) / (ss_w / df_w), ss_a / ss_t

    obs_f, obs_r2 = _pseudo_f(groups_arr)
    if np.isnan(obs_f):
        return np.nan, np.nan, np.nan

    perm_fs = []
    for _ in range(n_perm):
        pf, _ = _pseudo_f(np.random.permutation(groups_arr))
        if not np.isnan(pf):
            perm_fs.append(pf)

    p_val = (np.sum(np.array(perm_fs) >= obs_f) + 1) / (n_perm + 1) if perm_fs else np.nan

    print(f"    [*] PERMANOVA  F = {obs_f:.4f}   R² = {obs_r2:.4f}   p = {p_val:.4f}   (n_perm = {n_perm})")
    return float(obs_f), float(obs_r2), float(p_val)


def autopct_generator(pct):
    return f'{pct:.1f}%' if pct >= 1.5 else ''


def generate_global_pie_chart(df_rank, rank_level, threshold, output_base, fmt,
                               no_table=False):
    """Phase 1: Generates the Relative Abundance Probability Simplex."""
    print(f"[*] Extracting scalar probability simplex -> {output_base}")

    sample_cols = [col for col in df_rank.columns
                   if col not in ['Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name']]
    df_counts   = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

    rel_abund      = df_counts.div(df_counts.sum(axis=0), axis=1)
    mean_rel_abund = rel_abund.mean(axis=1)

    high_abund = mean_rel_abund[mean_rel_abund >= threshold]
    low_abund  = mean_rel_abund[mean_rel_abund < threshold]

    plot_data = high_abund.copy()
    if not low_abund.empty and low_abund.sum() > 0:
        plot_data['Others (<{:.1%})'.format(threshold)] = low_abund.sum()

    plot_data = plot_data.sort_values(ascending=False)

    fig, ax = plt.subplots(figsize=(14, 8))
    wedges, texts, autotexts = ax.pie(
        plot_data, labels=None, autopct=autopct_generator,
        startangle=140, colors=COLORS[:len(plot_data)],
        wedgeprops=dict(edgecolor='white', linewidth=1.5)
    )

    plt.setp(autotexts, size=11, weight="bold", color="white")
    for autotext in autotexts:
        autotext.set_path_effects([path_effects.withStroke(linewidth=2, foreground='black')])

    total_sum     = plot_data.sum()
    legend_labels = [f"{idx} ({ (val/total_sum)*100:.2f}% )"
                     for idx, val in plot_data.items()]

    ax.legend(wedges, legend_labels, title=f"{rank_level.capitalize()} Taxonomy",
              title_fontproperties={'weight': 'bold', 'size': 14}, loc="center left",
              bbox_to_anchor=(1.05, 0.5), fontsize=12, frameon=False,
              handlelength=1.5, handleheight=1.5)

    ax.set_title(f"Global Relative Abundance ({rank_level.capitalize()})",
                 fontsize=18, fontweight='bold', pad=20)
    plt.tight_layout()
    export_topological_projection(fig, output_base, fmt, pad_inches=0.5)
    plt.close(fig)

    # Pie chart summary table — Taxon | Mean_Rel_Abundance_pct
    if not no_table:
        table_df = pd.DataFrame({
            'Taxon':                  plot_data.index,
            'Mean_Rel_Abundance_pct': (plot_data.values / total_sum * 100).round(4)
        })
        export_summary_table(table_df, output_base, sheet_name='PieChart_Abundance')


def confidence_ellipse(x, y, ax, n_std=2.447, facecolor='none', **kwargs):
    """
    Mathematical Implementation of the Covariance Confidence Ellipse.
    n_std = 2.447 represents ~95% confidence for a 2D distribution.
    """
    if x.size < 3:
        return

    cov     = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])

    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)

    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x  = np.mean(x)
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y  = np.mean(y)

    transf = (transforms.Affine2D()
              .rotate_deg(45)
              .scale(scale_x, scale_y)
              .translate(mean_x, mean_y))

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def generate_pcoa_plot(df_rank, rank_level, metadata_path, category_col, sample_id_col,
                        output_base, fmt, unknown_mode='drop_all', no_table=False):
    """
    Phase 2: Generates the Principal Coordinate Analysis (PCoA) via Spectral Decomposition.
    Computes ANOSIM (Clarke 1993) and PERMANOVA (Anderson 2001) with 999 permutations
    and renders both R/F statistics + p-values in an annotation box (upper-left corner).

    unknown_mode controls how samples absent from metadata are handled:
      'drop_all'  — exclude from distance matrix AND plot (cleanest ordination)
      'drop_plot' — keep in distance matrix math, but omit from the final plot
      'keep'      — include in plot as an explicit 'Unknown' group
    """
    print(f"[*] Solving spectral decomposition for PCoA -> {output_base}")
    print(f"    [*] Unknown sample mode: '{unknown_mode}'")

    sample_cols = [col for col in df_rank.columns
                   if col not in ['Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name']]
    df_counts = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

    rel_abund    = df_counts.div(df_counts.sum(axis=0), axis=1).T
    sample_names = rel_abund.index.tolist()

    # --- Resolve metadata groups before touching the matrix ---
    groups = ["All Samples"] * len(sample_names)
    if metadata_path and category_col:
        meta_df = (pd.read_csv(metadata_path) if metadata_path.endswith('.csv')
                   else pd.read_excel(metadata_path))

        if pd.api.types.is_numeric_dtype(meta_df[category_col]):
            print(f"    [*] Continuous numerical variable detected for '{category_col}'.")
            print(f"    [*] Applying Quantile Binning (q=4) to ensure sufficient group sizes.")
            try:
                bins = pd.qcut(meta_df[category_col].dropna(), q=4)
                meta_df[category_col] = bins.astype(str)
            except ValueError:
                bins = pd.cut(meta_df[category_col].dropna(), bins=4)
                meta_df[category_col] = bins.astype(str)

        meta_df[sample_id_col] = meta_df[sample_id_col].astype(str).str.strip().str.lower()
        meta_map = dict(zip(meta_df[sample_id_col], meta_df[category_col].astype(str)))

        groups = []
        for name in sample_names:
            clean_name = name.split('_')[0].strip().lower()
            mapped_val = meta_map.get(clean_name, "Unknown")
            if mapped_val == 'nan':
                mapped_val = "Unknown"
            if mapped_val == "Unknown":
                print(f"    [!] Topological Mismatch: '{name}' (cleaned: '{clean_name}') NOT FOUND in metadata.")
            groups.append(mapped_val)

    # --- Apply unknown_mode: drop_all removes samples before matrix computation ---
    if unknown_mode == 'drop_all':
        keep_mask = [g != 'Unknown' for g in groups]
        n_dropped = keep_mask.count(False)
        if n_dropped > 0:
            dropped = [sample_names[i] for i, k in enumerate(keep_mask) if not k]
            print(f"    [!] drop_all: removing {n_dropped} unmatched sample(s): {dropped}")
        rel_abund    = rel_abund.iloc[[i for i, k in enumerate(keep_mask) if k]]
        sample_names = [s for s, k in zip(sample_names, keep_mask) if k]
        groups       = [g for g, k in zip(groups,       keep_mask) if k]
        if rel_abund.empty:
            print("    [!] No matched samples remaining after drop_all. Aborting PCoA plot.")
            return

    # Non-Euclidean Bray-Curtis metric calculation
    dist_array  = pdist(rel_abund.values, metric='braycurtis')
    dist_matrix = squareform(dist_array)

    # Torgerson-Gower scaling
    n = dist_matrix.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J.dot(dist_matrix ** 2).dot(J)

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    idx          = np.argsort(eigenvalues)[::-1]
    eigenvalues  = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    eigenvalues[eigenvalues < 0] = 0

    coords             = eigenvectors[:, :2] * np.sqrt(eigenvalues[:2])
    variance_explained = (eigenvalues / np.sum(eigenvalues)) * 100
    pc1_var, pc2_var   = variance_explained[0], variance_explained[1]

    df_pcoa = pd.DataFrame({
        'PC1': coords[:, 0], 'PC2': coords[:, 1],
        'Group': groups, 'Sample': sample_names
    })

    # --- Apply unknown_mode: drop_plot removes from visualization only ---
    if unknown_mode == 'drop_plot':
        unknown_mask = df_pcoa['Group'] == 'Unknown'
        n_dropped    = unknown_mask.sum()
        if n_dropped > 0:
            dropped = df_pcoa.loc[unknown_mask, 'Sample'].tolist()
            print(f"    [!] drop_plot: {n_dropped} sample(s) kept in matrix, hidden from plot: {dropped}")
        df_pcoa = df_pcoa[~unknown_mask].reset_index(drop=True)
        if df_pcoa.empty:
            print("    [!] No matched samples remaining after drop_plot. Aborting PCoA plot.")
            return

    # --- 'keep' mode: unknowns pass through as-is and are plotted as 'Unknown' group ---

    # --- Statistical Analysis: ANOSIM + PERMANOVA ---
    # Extract distance sub-matrix aligned with the final (post-filter) sample set.
    # This correctly handles all three unknown_mode strategies.
    _stats_samples  = df_pcoa['Sample'].tolist()
    _stats_groups   = df_pcoa['Group'].tolist()
    _sample_to_idx  = {s: i for i, s in enumerate(sample_names)}
    _pcoa_idx       = [_sample_to_idx[s] for s in _stats_samples if s in _sample_to_idx]
    _dist_stats     = dist_matrix[np.ix_(_pcoa_idx, _pcoa_idx)]

    anosim_r_val, anosim_p_val          = compute_anosim(_dist_stats, _stats_groups)
    perm_f_val,   perm_r2_val, perm_p_val = compute_permanova(_dist_stats, _stats_groups)

    # --- Plot ---
    unique_groups = df_pcoa['Group'].unique()

    fig, ax = plt.subplots(figsize=(11, 8), facecolor='#f8f9fa')
    ax.set_facecolor('#f4f4f6')
    ax.grid(color='white', linestyle='-', linewidth=1.5, alpha=0.8)

    for i, group in enumerate(unique_groups):
        subset = df_pcoa[df_pcoa['Group'] == group]
        color  = COLORS[i % len(COLORS)]

        ax.scatter(subset['PC1'], subset['PC2'],
                   s=180, alpha=0.9, label=group,
                   color=color, edgecolors='white', linewidth=2, zorder=3)

        if len(subset) >= 3 and group != "Unknown":
            confidence_ellipse(subset['PC1'].values, subset['PC2'].values, ax,
                               n_std=2.447, edgecolor=color, facecolor=color,
                               alpha=0.15, linewidth=2, zorder=2)
            confidence_ellipse(subset['PC1'].values, subset['PC2'].values, ax,
                               n_std=2.447, edgecolor=color, facecolor='none',
                               linestyle='--', alpha=0.8, linewidth=1.5, zorder=2)

    ax.axhline(0, color='black', linestyle=':', linewidth=1.2, alpha=0.6, zorder=1)
    ax.axvline(0, color='black', linestyle=':', linewidth=1.2, alpha=0.6, zorder=1)

    ax.set_xlabel(f"PC1 ({pc1_var:.1f}%)", fontsize=14, fontweight='bold', color='#333333')
    ax.set_ylabel(f"PC2 ({pc2_var:.1f}%)", fontsize=14, fontweight='bold', color='#333333')

    for spine in ax.spines.values():
        spine.set_visible(False)

    leg = ax.legend(title=category_col if category_col else "Group",
                    fontsize=12, title_fontproperties={'weight': 'bold', 'size': 13},
                    bbox_to_anchor=(1.03, 0.5), loc='center left', frameon=True,
                    facecolor='white', edgecolor='white', shadow=True)
    leg.get_title().set_color('#333333')

    # --- ANOSIM + PERMANOVA annotation box (upper-left corner) ---
    def _fmt(v, d=3):
        return f"{v:.{d}f}" if (v is not None and not np.isnan(v)) else "N/A"

    stats_text = (
        f"ANOSIM     R = {_fmt(anosim_r_val)}   p = {_fmt(anosim_p_val)}\n"
        f"PERMANOVA  F = {_fmt(perm_f_val, 2)}   R² = {_fmt(perm_r2_val)}   p = {_fmt(perm_p_val)}"
    )
    ax.text(0.03, 0.97, stats_text,
            transform=ax.transAxes, ha='left', va='top',
            fontsize=10, fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.92, linewidth=1.2),
            zorder=10)

    plt.tight_layout()
    export_topological_projection(fig, output_base, fmt, pad_inches=0.5)
    plt.close(fig)

    # PCoA summary table — Sheet 1: coordinates + groups; Sheet 2: Bray-Curtis matrix
    if not no_table:
        coords_df = df_pcoa[['Sample', 'Group', 'PC1', 'PC2']].round(6).reset_index(drop=True)

        dist_df = pd.DataFrame(
            _dist_stats,
            index=_stats_samples,
            columns=_stats_samples
        ).round(6)
        dist_df.insert(0, 'Sample', _stats_samples)

        export_summary_table(
            coords_df, output_base,
            sheet_name='PCoA_Coordinates',
            extra_sheets={'BrayCurtis_Distance': dist_df}
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PCoA and Pie Chart Topological Engine")
    parser.add_argument('-d',    '--data',      type=str, required=True,
                        help='Path to taxa Excel file.')
    parser.add_argument('-r',    '--rank',      type=str, default='genus',
                        help='Taxonomic rank.')
    parser.add_argument('-t',    '--threshold', type=float, default=0.01,
                        help='Abundance threshold for pie chart.')
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
    parser.add_argument('--mode', type=str, choices=['pie', 'pcoa', 'both'], default='both',
                        help='Routing logic: pie (simplex), pcoa (ordination), or both.')
    parser.add_argument('--unknown', type=str,
                        choices=['drop_all', 'drop_plot', 'keep'], default='drop_all',
                        help=(
                            'How to handle samples absent from metadata. '
                            'drop_all: remove from distance matrix AND plot. '
                            'drop_plot: keep in Bray-Curtis matrix, hide from plot. '
                            'keep: plot as an explicit "Unknown" group.'
                        ))
    parser.add_argument('--no_table', action='store_true',
                        help='Skip exporting summary tables (.xlsx).')

    args = parser.parse_args()

    OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'results', 'PCoA_PieCharts')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    safe_org_name = args.organism.replace(" ", "_")
    base_name     = args.output if args.output else f"{safe_org_name}_{args.rank}_{args.threshold}"
    output_base   = os.path.join(OUTPUT_DIR, base_name)

    # Data extraction phase
    df      = pd.read_excel(args.data, sheet_name=0)
    df_rank = df[df['Rank'].str.lower() == args.rank.lower()].copy()
    if df_rank.empty:
        raise ValueError(f"CRITICAL FAULT: No data found for the taxonomic level: {args.rank}")

    # Bifurcated operational execution
    if args.mode in ['pie', 'both']:
        generate_global_pie_chart(
            df_rank, args.rank, args.threshold,
            f"{output_base}_PieChart", args.format.lower(),
            no_table=args.no_table
        )

    if args.mode in ['pcoa', 'both']:
        if args.metadata and args.category:
            generate_pcoa_plot(
                df_rank, args.rank, args.metadata, args.category, args.sample_id,
                f"{output_base}_PCoA", args.format.lower(),
                unknown_mode=args.unknown,
                no_table=args.no_table
            )
        else:
            print("[!] WARNING: PCoA requires valid metadata and categorical vectors. Aborting sub-routine.")
