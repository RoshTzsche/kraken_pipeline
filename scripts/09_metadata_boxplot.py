"""
08_generate_lefse.py
====================
LEfSe — Linear Discriminant Analysis Effect Size
Anderson 2001 / Segata et al. 2011 (PLoS Comput Biol)

Pipeline-compatible script for the Kraken2 taxonomic classification output.
Now features automated combinatorial pairwise execution.
"""

import argparse
import os
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats as sp_stats
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from tqdm import tqdm

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Shared palette
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


def export_topological_projection(fig, output_base, fmt, pad_inches=0.4):
    if fmt == 'pdf':
        with PdfPages(f"{output_base}.pdf") as pdf:
            pdf.savefig(fig, bbox_inches='tight', pad_inches=pad_inches,
                        facecolor=fig.get_facecolor())
    else:
        ext = 'tif' if fmt == 'tiff' else 'png'
        fig.savefig(f"{output_base}.{ext}", format=fmt, dpi=300,
                    bbox_inches='tight', pad_inches=pad_inches,
                    facecolor=fig.get_facecolor())

def export_summary_table(df, output_base, sheet_name='LEfSe', extra_sheets=None):
    path = f"{output_base}_table.xlsx"
    with pd.ExcelWriter(path, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
        ws = writer.sheets[sheet_name]
        for col in ws.columns:
            w = max((len(str(c.value)) for c in col if c.value), default=8)
            ws.column_dimensions[col[0].column_letter].width = min(w + 4, 60)

        if extra_sheets:
            for sname, sdf in extra_sheets.items():
                sname = sname[:31]
                sdf.to_excel(writer, index=False, sheet_name=sname)
                ws2 = writer.sheets[sname]
                for col in ws2.columns:
                    w = max((len(str(c.value)) for c in col if c.value), default=8)
                    ws2.column_dimensions[col[0].column_letter].width = min(w + 4, 60)
    print(f"    [*] Table saved: {path}")

def relative_abundance(df_counts):
    totals = df_counts.sum(axis=0)
    totals = totals.replace(0, np.nan)
    return df_counts.div(totals, axis=1)

def kruskal_wallis_test(ra_df, groups, alpha=0.05):
    unique_groups = sorted(set(groups))
    passing = {}
    for feat in tqdm(ra_df.index, desc='  Kruskal-Wallis', unit='feat', ncols=72, colour='cyan'):
        grp_arrays = [ra_df.loc[feat, [s for s, g in zip(ra_df.columns, groups) if g == ug]].values
                      for ug in unique_groups]
        valid = [a for a in grp_arrays if len(a) > 0 and a.std() >= 0]
        if len(valid) < 2:
            continue
        try:
            H, p = sp_stats.kruskal(*valid)
        except ValueError:
            continue
        if p < alpha:
            passing[feat] = (H, p)
    return passing

def wilcoxon_pairwise_test(ra_df, groups, candidate_features, alpha=0.05):
    unique_groups = sorted(set(groups))
    groups_arr    = np.array(groups)
    passed = {}
    for feat in tqdm(candidate_features, desc='  Wilcoxon      ', unit='feat', ncols=72, colour='green'):
        row        = ra_df.loc[feat]
        medians    = {ug: row[[s for s, g in zip(ra_df.columns, groups) if g == ug]].median()
                      for ug in unique_groups}
        winner     = max(medians, key=medians.get)

        significant = False
        for i in range(len(unique_groups)):
            for j in range(i + 1, len(unique_groups)):
                g1, g2 = unique_groups[i], unique_groups[j]
                v1 = row[groups_arr == g1].values
                v2 = row[groups_arr == g2].values
                if len(v1) < 2 or len(v2) < 2:
                    continue
                try:
                    _, p = sp_stats.mannwhitneyu(v1, v2, alternative='two-sided')
                    if p < alpha:
                        significant = True
                        break
                except ValueError:
                    continue
            if significant:
                break
        if significant:
            passed[feat] = winner
    return passed

def compute_lda_scores(ra_df, groups, candidate_features, winning_classes):
    groups_arr    = np.array(groups)
    scores        = {}
    X_full = np.log1p(ra_df.loc[candidate_features].T.values)   

    for feat_idx, feat in enumerate(tqdm(candidate_features, desc='  LDA scoring   ', unit='feat', ncols=72, colour='yellow')):
        x = X_full[:, feat_idx].reshape(-1, 1)   
        winner = winning_classes[feat]
        y_bin = (groups_arr == winner).astype(int)

        if len(np.unique(y_bin)) < 2:
            scores[feat] = 0.0
            continue

        try:
            lda = LinearDiscriminantAnalysis(solver='svd')
            lda.fit(x, y_bin)
            ld1 = lda.transform(x).ravel()
            mean_winner = ld1[y_bin == 1].mean()
            mean_rest   = ld1[y_bin == 0].mean()
            raw_score   = np.log10(abs(mean_winner - mean_rest) + 1e-9)
            raw_score = max(abs(raw_score), 0)
            scores[feat] = raw_score if mean_winner > mean_rest else -raw_score
        except Exception:
            scores[feat] = 0.0

    return scores

def run_lefse(df_rank, rank_level, metadata_path, category_col, sample_id_col,
              output_base, fmt, lda_threshold=2.0, kw_alpha=0.05, wilcox_alpha=0.05,
              top_n=3, sort_by_lda=False, label_col=None, no_table=False,
              focus_groups=None):

    comp_name = f" ({focus_groups[0]} vs {focus_groups[1]})" if focus_groups else ""
    print(f"\n[*] Running LEfSe{comp_name} → {output_base}")

    df_work = df_rank.copy()
    if label_col and label_col in df_work.columns:
        df_work.index = df_work[label_col].astype(str)
    elif 'Name' in df_work.columns:
        df_work.index = df_work['Name'].astype(str)
    elif 'Scientific Name' in df_work.columns:
        df_work.index = df_work['Scientific Name'].astype(str)
    else:
        df_work.index = df_work.index.astype(str)

    meta_cols = [c for c in ['Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name']
                 if c in df_work.columns]
    sample_cols = [c for c in df_work.columns if c not in meta_cols]

    df_counts = df_work[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

    if metadata_path.endswith('.csv'):
        meta = pd.read_csv(metadata_path)
    else:
        meta = pd.read_excel(metadata_path)

    meta[sample_id_col] = meta[sample_id_col].astype(str).str.strip().str.lower()
    meta_map = dict(zip(meta[sample_id_col], meta[category_col].astype(str)))

    all_samples   = list(df_counts.columns)
    meta_keys     = list(meta_map.keys())
    groups_raw    = []
    for name in all_samples:
        clean = name.split('_')[0].strip().lower()
        val   = meta_map.get(clean, 'Unknown')

        if val in ('Unknown', 'nan'):
            candidates = [k for k in meta_keys if clean.startswith(k)]
            if candidates:
                best_key = max(candidates, key=len)
                val      = meta_map[best_key]

        if val == 'nan':
            val = 'Unknown'
        groups_raw.append(val)

    # DYNAMIC FILTERING: If focus_groups are provided, strictly limit the analysis to those groups
    if focus_groups:
        keep_mask = [g in focus_groups for g in groups_raw]
    else:
        keep_mask = [g != 'Unknown' for g in groups_raw]

    common_samples = [s for s, k in zip(all_samples, keep_mask) if k]
    groups         = [g for g, k in zip(groups_raw,  keep_mask) if k]

    if len(common_samples) < 3:
        print(f"    [!] CRITICAL: Not enough samples matched for {comp_name}. Skipping.")
        return

    df_counts = df_counts[common_samples]
    unique_groups = sorted(set(groups))

    ra_df = relative_abundance(df_counts).fillna(0)
    ra_df = ra_df.loc[ra_df.var(axis=1) > 0]

    kw_results = kruskal_wallis_test(ra_df, groups, alpha=kw_alpha)
    if not kw_results:
        print("    [!] No features passed Kruskal-Wallis. Exiting.")
        return

    winning_classes = wilcoxon_pairwise_test(ra_df, groups, list(kw_results.keys()), alpha=wilcox_alpha)
    if not winning_classes:
        print("    [!] No features passed Wilcoxon. Exiting.")
        return

    lda_scores = compute_lda_scores(ra_df, groups, list(winning_classes.keys()), winning_classes)

    records = []
    for feat, winner in winning_classes.items():
        H, kw_p = kw_results[feat]
        lda      = lda_scores.get(feat, 0.0)
        abs_lda  = abs(lda)
        if abs_lda < lda_threshold:
            continue
        records.append({
            'Feature':       feat,
            'Class':         winner,
            'LDA_score':     round(lda,     4),
            'Abs_LDA':       round(abs_lda, 4),
            'KW_H':          round(H,       4),
            'KW_p':          round(kw_p,    6),
            'Mean_RA':       round(float(ra_df.loc[feat].mean()), 6),
        })

    if not records:
        print(f"    [!] No features passed |LDA| >= {lda_threshold}. Try lowering --lda_threshold.")
        return

    df_results = pd.DataFrame(records).sort_values('Abs_LDA', ascending=False)
    
    if top_n > 0:
        df_plot = df_results.sort_values('Abs_LDA', ascending=False)
        df_plot = df_plot.groupby('Class').head(top_n).copy()
    else:
        df_plot = df_results.copy()

    group_means = []
    for feat in df_plot['Feature']:
        row = {'Feature': feat}
        for ug in unique_groups:
            samp = [s for s, g in zip(common_samples, groups) if g == ug]
            row[f'MeanRA_{ug}'] = round(float(ra_df.loc[feat, samp].mean()), 6)
        group_means.append(row)
    df_group_means = pd.DataFrame(group_means)

    # Assign distinct colors per unique group for consistency
    color_map = {g: COLORS[i % len(COLORS)] for i, g in enumerate(unique_groups)}

    if sort_by_lda:
        df_plot = df_plot.sort_values('LDA_score')
    else:
        df_plot = df_plot.sort_values(['Class', 'Abs_LDA'], ascending=[True, False])

    n_bars   = len(df_plot)
    fig_h    = max(5, n_bars * 0.38 + 2)
    fig, ax  = plt.subplots(figsize=(11, fig_h))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('#F7F7F7')

    for i, row in enumerate(df_plot.itertuples()):
        lda_val = row.LDA_score
        color   = color_map[row.Class]
        bar_val = lda_val if sort_by_lda else row.Abs_LDA
        bar_dir = 1 if lda_val >= 0 else -1

        if not sort_by_lda:
            ax.barh(i, bar_val, color=color, edgecolor='white', linewidth=0.5, height=0.72, zorder=3)
        else:
            ax.barh(i, lda_val, color=color, edgecolor='white', linewidth=0.5, height=0.72, zorder=3)

        label_x = bar_val + 0.05 if not sort_by_lda else (lda_val + 0.05 * bar_dir)
        ax.text(label_x, i, f"{row.Abs_LDA:.2f}", va='center', ha='left', fontsize=8.5, color='#333333', zorder=4)

    ax.set_yticks(range(n_bars))
    ax.set_yticklabels(df_plot['Feature'].tolist(), fontsize=9.5)

    ax.axvline(0, color='#888888', linewidth=0.8, linestyle='--', zorder=2)
    ax.set_xlabel('LDA Score (log₁₀)', fontsize=12, labelpad=8)
    ax.set_title(f'LEfSe — {rank_level.capitalize()} level{comp_name}\n[|LDA| ≥ {lda_threshold} | KW p<{kw_alpha} | Wilcoxon p<{wilcox_alpha}]',
                 fontsize=13, fontweight='bold', pad=14)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#CCCCCC')
    ax.tick_params(axis='both', which='both', length=0)
    ax.grid(axis='x', color='white', linewidth=1.2, zorder=1)
    ax.set_xlim(left=0 if not sort_by_lda else None)

    patches = [mpatches.Patch(color=color_map[g], label=g) for g in unique_groups]
    ax.legend(handles=patches, title='Class', title_fontsize=10, fontsize=9.5, loc='lower right',
              frameon=True, framealpha=0.9, edgecolor='#CCCCCC')

    subtitle_text = f"n = {len(df_results)} total significant biomarkers"
    if top_n > 0:
        subtitle_text += f" | showing top {top_n} per class"
    else:
        subtitle_text += f" | showing all"
        
    ax.text(0.01, 1.01, subtitle_text,
            transform=ax.transAxes, fontsize=9, color='#666666', ha='left')

    plt.tight_layout(pad=1.5)
    export_topological_projection(fig, output_base, fmt)
    plt.close(fig)

    if not no_table:
        export_summary_table(df_results, output_base, sheet_name='LEfSe_Biomarkers',
                             extra_sheets={'Group_Mean_RelAbund': df_group_means})


def main():
    parser = argparse.ArgumentParser(description='LEfSe — Linear Discriminant Analysis Effect Size.')
    parser.add_argument('-d',  '--data',          type=str, required=True)
    parser.add_argument('-r',  '--rank',           type=str, required=True)
    parser.add_argument('-m',  '--metadata',       type=str, required=True)
    parser.add_argument('-c',  '--category',       type=str, required=True)
    parser.add_argument('-id', '--sample_id',      type=str, default='SampleID')
    parser.add_argument('-fmt','--format',         type=str, choices=['pdf', 'png', 'tiff'], default='pdf')
    parser.add_argument('-o',  '--output',         type=str, default=None)
    parser.add_argument('-org','--organism',       type=str, default='Microbiome')
    parser.add_argument('--lda_threshold',         type=float, default=2.0)
    parser.add_argument('--kw_alpha',              type=float, default=0.05)
    parser.add_argument('--wilcox_alpha',          type=float, default=0.05)
    parser.add_argument('--top',                   type=int,   default=3,
                        help='Maximum number of biomarkers to display PER CATEGORY (default: 3. Set to 0 for ALL).')
    parser.add_argument('--rank_by_lda',           action='store_true')
    parser.add_argument('--label_col',             type=str, default=None)
    parser.add_argument('--no_table',              action='store_true')
    
    # NEW ARGUMENT
    parser.add_argument('--pairwise',              action='store_true',
                        help='Run exhaustive pairwise combinations between all classes in the metadata category.')

    args = parser.parse_args()

    OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'LEfSe')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    safe_org  = args.organism.replace(' ', '_')
    base_name = args.output if args.output else f"{safe_org}_{args.rank}_LEfSe"
    
    df      = pd.read_excel(args.data, sheet_name=0)
    df_rank = df[df['Rank'].str.lower() == args.rank.lower()].copy()

    if df_rank.empty:
        raise ValueError(f"No data found for rank: {args.rank}")

    # Determine execution mode
    if args.pairwise:
        # Pre-read metadata to find unique groups
        if args.metadata.endswith('.csv'):
            meta = pd.read_csv(args.metadata)
        else:
            meta = pd.read_excel(args.metadata)
            
        unique_classes = meta[args.category].dropna().astype(str).unique().tolist()
        unique_classes = [c for c in unique_classes if c.lower() != 'nan']
        
        if len(unique_classes) < 2:
            raise ValueError("Pairwise mode requires at least 2 distinct classes in the metadata.")
            
        combos = list(combinations(sorted(unique_classes), 2))
        print(f"[*] Pairwise Mode Enabled: Found {len(unique_classes)} classes, running {len(combos)} combinations.")
        
        for g1, g2 in combos:
            # Modify output base to avoid overwriting files
            pair_base = os.path.join(OUTPUT_DIR, f"{base_name}_{g1}_vs_{g2}")
            
            run_lefse(
                df_rank       = df_rank,
                rank_level    = args.rank,
                metadata_path = args.metadata,
                category_col  = args.category,
                sample_id_col = args.sample_id,
                output_base   = pair_base,
                fmt           = args.format.lower(),
                lda_threshold = args.lda_threshold,
                kw_alpha      = args.kw_alpha,
                wilcox_alpha  = args.wilcox_alpha,
                top_n         = args.top,
                sort_by_lda   = args.rank_by_lda,
                label_col     = args.label_col,
                no_table      = args.no_table,
                focus_groups  = [g1, g2] # Pass the constrained groups
            )
            
    else:
        # Standard Multi-Class Mode
        output_base = os.path.join(OUTPUT_DIR, base_name)
        run_lefse(
            df_rank       = df_rank,
            rank_level    = args.rank,
            metadata_path = args.metadata,
            category_col  = args.category,
            sample_id_col = args.sample_id,
            output_base   = output_base,
            fmt           = args.format.lower(),
            lda_threshold = args.lda_threshold,
            kw_alpha      = args.kw_alpha,
            wilcox_alpha  = args.wilcox_alpha,
            top_n         = args.top,
            sort_by_lda   = args.rank_by_lda,
            label_col     = args.label_col,
            no_table      = args.no_table,
            focus_groups  = None
        )

if __name__ == '__main__':
    main()
