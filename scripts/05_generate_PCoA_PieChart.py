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
from statsmodels.stats.multitest import multipletests
from analysis_support import (count_matrix, metadata_groups, group_order, read_frame,
                              write_run_record, pcoa_lingoes, dispersion_test)

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
    from skbio import DistanceMatrix
    from skbio.stats.distance import anosim
    if len(set(groups)) < 2 or min(pd.Series(groups).value_counts()) < 2:
        return np.nan, np.nan
    result = anosim(DistanceMatrix(dist_matrix), list(groups), permutations=n_perm, seed=seed)
    return float(result['test statistic']), float(result['p-value'])


def compute_permanova(dist_matrix, groups, n_perm=999, seed=42):
    from skbio import DistanceMatrix
    from skbio.stats.distance import permanova
    if len(set(groups)) < 2 or min(pd.Series(groups).value_counts()) < 2:
        return np.nan, np.nan, np.nan
    result = permanova(DistanceMatrix(dist_matrix), list(groups), permutations=n_perm, seed=seed)
    f = float(result['test statistic'])
    k, n = len(set(groups)), len(groups)
    r2 = f * (k - 1) / (f * (k - 1) + n - k) if np.isfinite(f) else (1.0 if np.isinf(f) else np.nan)
    return f, r2, float(result['p-value'])


def compute_pairwise_anosim(dist_matrix, groups, n_perm=999, seed=42):
    from itertools import combinations
    groups = np.asarray(groups)
    rows = []
    for a, b in combinations(group_order(groups), 2):
        idx = np.flatnonzero((groups == a) | (groups == b))
        r, p = compute_anosim(dist_matrix[np.ix_(idx, idx)], groups[idx], n_perm, seed)
        rows.append({'Group1': a, 'Group2': b, 'R': r, 'p_value': p})
    frame = pd.DataFrame(rows, columns=['Group1', 'Group2', 'R', 'p_value'])
    frame['p_adjusted'] = np.nan
    finite = frame.p_value.notna()
    if finite.any():
        frame.loc[finite, 'p_adjusted'] = multipletests(frame.loc[finite, 'p_value'], method='fdr_bh')[1]
    frame['Adjustment'] = 'BH across pairwise ANOSIM comparisons'
    return frame


def autopct_generator(pct):
    # Left for compatibility, not used in the new 3D pie chart callouts
    return f'{pct:.1f}%' if pct >= 1.5 else ''


def generate_global_pie_chart(df_rank, rank_level, threshold, output_base, fmt,
                               no_table=False, label_col=None):
    """Phase 1: Generates the Relative Abundance Probability Simplex with 3D Callouts."""
    print(f"[*] Extracting scalar probability simplex -> {output_base}")

    sample_cols = [col for col in df_rank.columns
                   if col not in ['Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name']]

    # --- Resolve label index: use label_col if provided and present, else fall back to TaxID/Name ---
    df_work = df_rank.copy()
    if label_col and label_col in df_work.columns:
        print(f"    [*] Pie chart labels resolved from column: '{label_col}'")
        df_work.index = df_work[label_col].astype(str)
    elif 'Name' in df_work.columns:
        df_work.index = df_work['Name'].astype(str)
    elif 'Scientific Name' in df_work.columns:
        df_work.index = df_work['Scientific Name'].astype(str)
    else:
        print("    [!] No label column found — using row index as label.")

    df_counts = count_matrix(df_work).groupby(level=0).sum()

    rel_abund      = df_counts.div(df_counts.sum(axis=0).replace(0, np.nan), axis=1)
    mean_rel_abund = rel_abund.mean(axis=1)

    high_abund = mean_rel_abund[mean_rel_abund >= threshold]
    low_abund  = mean_rel_abund[mean_rel_abund < threshold]

    plot_data = high_abund.copy()
    if not low_abund.empty and low_abund.sum() > 0:
        plot_data['Others (<{:.1%})'.format(threshold)] = low_abund.sum()

    plot_data = plot_data.sort_values(ascending=False)
    total_sum = plot_data.sum()

    fig, ax = plt.subplots(figsize=(16, 9))
    
    explode = [0.03] * len(plot_data)

    wedges, texts = ax.pie(
        plot_data, labels=None, startangle=140, 
        colors=COLORS[:len(plot_data)], explode=explode, shadow=False,
        wedgeprops=dict(edgecolor='white', linewidth=1.0)
    )

    LABEL_R   = 1.15
    TEXT_X    = 1.55
    MIN_SEP   = 0.13
    FONT_SIZE = 11

    label_info = []
    for i, wedge in enumerate(wedges):
        ang      = (wedge.theta2 - wedge.theta1) / 2.0 + wedge.theta1
        ang_rad  = np.deg2rad(ang)
        xi       = np.cos(ang_rad)
        yi       = np.sin(ang_rad)
        pct      = (plot_data.values[i] / total_sum) * 100
        label_info.append({
            'wedge':   wedge,
            'ang':     ang,
            'xi':      xi,
            'yi':      yi,
            'bend_x':  LABEL_R * xi,
            'bend_y':  LABEL_R * yi,
            'side':    'right' if xi >= 0 else 'left',
            'text':    f"{plot_data.index[i]} ({pct:.1f}%)",
        })

    def _resolve_collisions(items, min_sep):
        if not items:
            return items
        items = sorted(items, key=lambda d: d['bend_y'])
        if len(items) == 1:
            items[0]['label_y'] = items[0]['bend_y']
            return items
        placed = [d['bend_y'] for d in items]

        for _ in range(50):
            moved = False
            for k in range(1, len(placed)):
                gap = placed[k] - placed[k - 1]
                if gap < min_sep:
                    shift = (min_sep - gap) / 2.0
                    placed[k - 1] -= shift
                    placed[k]     += shift
                    moved = True
            if not moved:
                break

        for d, y_new in zip(items, placed):
            d['label_y'] = y_new
        return items

    right_items = _resolve_collisions([d for d in label_info if d['side'] == 'right'], MIN_SEP)
    left_items  = _resolve_collisions([d for d in label_info if d['side'] == 'left'],  MIN_SEP)

    for d in right_items + left_items:
        ha      = 'left' if d['side'] == 'right' else 'right'
        text_x  = TEXT_X  if d['side'] == 'right' else -TEXT_X
        label_y = d['label_y']

        ax.plot(
            [d['xi'],     d['bend_x'], text_x],
            [d['yi'],     d['bend_y'], label_y],
            color='#444444', linewidth=0.9, zorder=4, solid_capstyle='round'
        )

        ax.text(
            text_x, label_y, d['text'],
            ha=ha, va='center',
            fontsize=FONT_SIZE, fontweight='bold',
            zorder=5,
        )

    ax.set_title(f"Global Relative Abundance ({rank_level.capitalize()})",
                 fontsize=18, fontweight='bold', pad=30)
                 
    plt.tight_layout()
    export_topological_projection(fig, output_base, fmt, pad_inches=0.5)
    plt.close(fig)

    if not no_table:
        table_df = pd.DataFrame({
            'Taxon':                  plot_data.index,
            'Mean_Rel_Abundance_pct': (plot_data.values / total_sum * 100).round(4)
        })
        export_summary_table(table_df, output_base, sheet_name='PieChart_Abundance')


def confidence_ellipse(x, y, ax, n_std=2.447, facecolor='none', **kwargs):
    """
    Approximate 95% Gaussian covariance ellipse for the distribution of points.
    It is a descriptive spread ellipse, not a confidence region for the mean.
    """
    if x.size < 3:
        return

    cov     = np.cov(x, y)
    if not np.isfinite(cov).all() or cov[0, 0] <= 0 or cov[1, 1] <= 0:
        return
    pearson = np.clip(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]), -1, 1)

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
                       output_base, fmt, unknown_mode='drop_all', no_table=False,
                       permutations=999, seed=42, data_path=None):
    """Bray-Curtis inference and Lingoes-corrected PCoA; unmatched groups never tested."""
    if permutations < 1:
        raise ValueError('At least one permutation is required')
    counts = count_matrix(df_rank)
    mapping, audit = metadata_groups(counts.columns, metadata_path, category_col, sample_id_col)
    totals = counts.sum(axis=0)
    empty = audit.Sample.isin(totals[totals.eq(0)].index)
    audit.loc[empty, ['Included', 'Reason']] = [False, 'no_taxonomic_detections']
    samples = [s for s in counts if totals[s] > 0 and (unknown_mode != 'drop_all' or mapping[s] is not None)]
    if len(samples) < 3:
        raise ValueError('PCoA needs at least three nonempty samples after metadata filtering')
    relative = counts[samples].div(totals[samples], axis=1).T
    distances = squareform(pdist(relative.to_numpy(), metric='braycurtis'))
    coords, explained, _, diagnostic = pcoa_lingoes(distances)
    coordinates = pd.DataFrame({'Sample': samples, 'Group': [mapping[s] or 'Unknown' for s in samples],
                                'PCo1': coords[:, 0], 'PCo2': coords[:, 1]})
    known = [i for i, sample in enumerate(samples) if mapping[sample] is not None]
    test_samples = [samples[i] for i in known]
    groups = [mapping[s] for s in test_samples]
    test_distances = distances[np.ix_(known, known)]
    valid = len(set(groups)) >= 2 and min(pd.Series(groups).value_counts()) >= 2
    if valid and np.any(test_distances):
        ar, ap = compute_anosim(test_distances, groups, permutations, seed)
        f, r2, pp = compute_permanova(test_distances, groups, permutations, seed)
        dispersion = dispersion_test(test_distances, test_samples, groups, permutations, seed)
        status = 'computed'
    else:
        ar = ap = f = r2 = pp = np.nan
        status = 'insufficient_replication' if not valid else 'all_distances_zero'
        dispersion = {'Test': 'PERMDISP', 'Statistic': np.nan, 'p_value': np.nan,
                      'Center': 'spatial median', 'Distance': 'Lingoes-corrected Bray-Curtis'}
    global_tests = pd.DataFrame([
        {'Test': 'PERMANOVA', 'Statistic': f, 'R_squared': r2, 'p_value': pp, 'Distance': 'Bray-Curtis'},
        {'Test': 'ANOSIM (secondary)', 'Statistic': ar, 'p_value': ap, 'Distance': 'Bray-Curtis'},
        dispersion,
    ])
    global_tests['Status'] = status
    global_tests['Permutations'] = permutations
    global_tests['Seed'] = seed
    global_tests['N_samples'] = len(groups)
    global_tests['N_groups'] = len(set(groups))
    pairwise = compute_pairwise_anosim(test_distances, groups, permutations, seed) if valid and np.any(test_distances) else pd.DataFrame()
    distance_frame = pd.DataFrame(distances, index=samples, columns=samples)
    coordinates.to_csv(f'{output_base}_coordinates.csv', index=False)
    distance_frame.to_csv(f'{output_base}_distances.csv', index_label='Sample')
    global_tests.to_csv(f'{output_base}_global_tests.csv', index=False)
    pairwise.to_csv(f'{output_base}_pairwise_ANOSIM.csv', index=False)
    write_run_record(output_base, 'beta_diversity', inputs=[data_path, metadata_path], audit=audit,
        parameters={'rank': rank_level, 'category': category_col, 'permutations': permutations,
                    'seed': seed, 'permutations_scheme': 'unrestricted; independent specimens required',
                    'distance': 'Bray-Curtis on per-sample proportions', 'ordination': 'PCoA with Lingoes correction',
                    'dispersion_center': 'spatial median', 'unknown_mode': unknown_mode,
                    'unknown_inference': 'always excluded', **diagnostic})
    plot_data = coordinates if unknown_mode == 'keep' else coordinates[coordinates.Group.ne('Unknown')]
    fig, ax = plt.subplots(figsize=(10, 7), facecolor='white')
    for i, group in enumerate(group_order(plot_data.Group)):
        subset = plot_data[plot_data.Group.eq(group)]
        color = COLORS[i % len(COLORS)]
        ax.scatter(subset.PCo1, subset.PCo2, s=70, label=f'{group} (n={len(subset)})',
                   color=color, edgecolors='black', linewidths=0.5, zorder=3)
        if len(subset) >= 3 and group != 'Unknown':
            confidence_ellipse(subset.PCo1.to_numpy(), subset.PCo2.to_numpy(), ax,
                               edgecolor=color, facecolor=color, alpha=0.12)
    ax.axhline(0, color='grey', linewidth=0.5)
    ax.axvline(0, color='grey', linewidth=0.5)
    ax.set_xlabel(f'PCo1 ({explained[0]:.1f}%)')
    ax.set_ylabel(f'PCo2 ({explained[1]:.1f}%)')
    ax.legend(title=category_col, bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
    ax.set_title(f'{rank_level.capitalize()} Bray-Curtis PCoA')
    if status == 'computed':
        ax.text(0.02, 0.98, f'PERMANOVA: R²={r2:.3f}, p={pp:.3g}\nPERMDISP: p={dispersion["p_value"]:.3g}',
                transform=ax.transAxes, va='top', fontsize=9)
    fig.tight_layout()
    export_topological_projection(fig, output_base, fmt)
    plt.close(fig)
    if not no_table:
        export_summary_table(coordinates, output_base, 'PCoA_Coordinates', {
            'BrayCurtis_Distance': distance_frame.reset_index(names='Sample'),
            'Global_Tests': global_tests, 'Pairwise_ANOSIM': pairwise})
    return coordinates, global_tests


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
    parser.add_argument('--label_col', type=str, default=None,
                        help='Column to use as pie chart slice labels (e.g. "Name", "Scientific Name"). '
                             'Defaults to "Name" if present, then "Scientific Name", then row index.')
    parser.add_argument('--no_table', action='store_true',
                        help='Skip exporting summary tables (.xlsx).')

    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'results', 'PCoA_PieCharts')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    safe_org_name = args.organism.replace(" ", "_")
    base_name     = args.output if args.output else f"{safe_org_name}_{args.rank}_{args.threshold}"
    output_base   = os.path.join(OUTPUT_DIR, base_name)

    # Data extraction phase
    df      = read_frame(args.data)
    df_rank = df[df['Rank'].str.lower() == args.rank.lower()].copy()
    if df_rank.empty:
        raise ValueError(f"CRITICAL FAULT: No data found for the taxonomic level: {args.rank}")

    # Bifurcated operational execution
    if args.mode in ['pie', 'both']:
        generate_global_pie_chart(
            df_rank, args.rank, args.threshold,
            f"{output_base}_PieChart", args.format.lower(),
            no_table=args.no_table,
            label_col=args.label_col
        )

    if args.mode in ['pcoa', 'both']:
        if args.metadata and args.category:
            generate_pcoa_plot(
                df_rank, args.rank, args.metadata, args.category, args.sample_id,
                f"{output_base}_PCoA", args.format.lower(),
                unknown_mode=args.unknown,
                no_table=args.no_table, permutations=args.permutations, seed=args.seed, data_path=args.data
            )
        else:
            print("[!] WARNING: PCoA requires valid metadata and categorical vectors. Aborting sub-routine.")

