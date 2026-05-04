"""
06_generate_LEfSe.py
====================
LEfSe — Linear Discriminant Analysis Effect Size
Anderson 2001 / Segata et al. 2011 (PLoS Comput Biol)

Pipeline-compatible script for the Kraken2 taxonomic classification output.

Usage
-----
python 06_generate_LEfSe.py \
    -d  ../results/final_tables/taxonomic_classification_clean.xlsx \
    -r  genus \
    -m  metadata.csv \
    -c  Month \
    -fmt png

Optional
--------
--rank_by_lda     Sort bars by LDA score (default: grouped by winning class)
--top             Number of top taxa to show (default: 30)
--lda_threshold   Minimum |LDA| score to include (default: 2.0)
--kw_alpha        Kruskal-Wallis significance cut-off (default: 0.05)
--wilcox_alpha    Wilcoxon significance cut-off (default: 0.05)
--no_table        Skip Excel export
--label_col       Column to use as feature label (default: "Name")
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats as sp_stats
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Shared palette — identical across pipeline scripts
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


# ===========================================================================
# I/O helpers
# ===========================================================================

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


# ===========================================================================
# Statistical core
# ===========================================================================

def relative_abundance(df_counts):
    """Row-normalise to relative abundance (proportions)."""
    totals = df_counts.sum(axis=0)
    totals = totals.replace(0, np.nan)
    return df_counts.div(totals, axis=1)


def kruskal_wallis_test(ra_df, groups, alpha=0.05):
    """
    Run Kruskal-Wallis H-test for each feature across groups.
    Returns a dict: feature → (H, p)  only for features passing alpha.
    """
    unique_groups = sorted(set(groups))
    passing = {}
    for feat in ra_df.index:
        grp_arrays = [ra_df.loc[feat, [s for s, g in zip(ra_df.columns, groups) if g == ug]].values
                      for ug in unique_groups]
        # need at least 2 non-degenerate groups
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
    """
    For each candidate feature, run pairwise Wilcoxon rank-sum tests
    between every group pair. A feature passes if AT LEAST ONE pair
    is significant (consistent with LEfSe's all-against-one strategy).

    Returns dict: feature → winning_class (group with highest median)
    """
    unique_groups = sorted(set(groups))
    groups_arr    = np.array(groups)
    passed = {}

    for feat in candidate_features:
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
    """
    Fit a one-vs-rest LDA per feature using the log-transformed relative
    abundance. Returns dict: feature → signed LDA score.

    Sign convention: positive = enriched in winning_class.
    """
    groups_arr    = np.array(groups)
    unique_groups = sorted(set(groups))
    scores        = {}

    # Build full feature matrix (samples × features) log-transformed
    X_full = np.log1p(ra_df.loc[candidate_features].T.values)   # (n_samples, n_feats)

    for feat_idx, feat in enumerate(candidate_features):
        x = X_full[:, feat_idx].reshape(-1, 1)   # single feature

        winner = winning_classes[feat]

        # Binary labels: winner vs rest
        y_bin = (groups_arr == winner).astype(int)

        if len(np.unique(y_bin)) < 2:
            scores[feat] = 0.0
            continue

        try:
            lda = LinearDiscriminantAnalysis(solver='svd')
            lda.fit(x, y_bin)
            # LDA score = absolute value of the mean difference in LD1 space
            ld1 = lda.transform(x).ravel()
            mean_winner = ld1[y_bin == 1].mean()
            mean_rest   = ld1[y_bin == 0].mean()
            raw_score   = np.log10(abs(mean_winner - mean_rest) + 1e-9)
            # Scale to typical LEfSe range (~2–5)
            raw_score = max(abs(raw_score), 0)
            scores[feat] = raw_score if mean_winner > mean_rest else -raw_score
        except Exception:
            scores[feat] = 0.0

    return scores


# ===========================================================================
# Cladogram
# ===========================================================================

RANK_ORDER = ['domain', 'superkingdom', 'phylum', 'class',
              'order', 'family', 'genus', 'species']

def _get_label(row, label_col):
    """Pick best display name from a DataFrame row."""
    if label_col and label_col in row.index and pd.notna(row[label_col]):
        return str(row[label_col]).strip()
    for col in ('Name', 'Scientific Name', 'TaxID'):
        if col in row.index and pd.notna(row[col]):
            return str(row[col]).strip()
    return str(row.name)


def build_taxonomy_tree(df_full, common_samples, label_col=None):
    """
    Build a lightweight tree from all ranks in the data.

    Returns
    -------
    nodes : dict  node_id -> {name, rank, rank_idx, children, parent,
                               counts (np.array over common_samples)}
    present_ranks : list of rank strings in RANK_ORDER order
    """
    meta_cols = {'Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name'}

    # Keep only columns that are actual samples
    available_samples = [c for c in common_samples if c in df_full.columns]

    # Determine which ranks are present
    all_ranks_lc = df_full['Rank'].str.lower().unique()
    present_ranks = [r for r in RANK_ORDER if r in all_ranks_lc]

    if not present_ranks:
        return {}, []

    nodes = {}
    # Virtual root
    nodes['__root__'] = dict(name='root', rank='root', rank_idx=-1,
                              children=[], parent=None,
                              counts=np.zeros(len(available_samples)))

    for r_idx, rank in enumerate(present_ranks):
        subset = df_full[df_full['Rank'].str.lower() == rank].copy()
        for _, row in subset.iterrows():
            label  = _get_label(row, label_col)
            tid    = str(row['TaxID']) if 'TaxID' in row.index else label
            nid    = f"{rank}::{tid}::{label}"
            counts = np.array([
                float(row[s]) if s in row.index else 0.0
                for s in available_samples
            ], dtype=float)
            nodes[nid] = dict(name=label, rank=rank, rank_idx=r_idx,
                               children=[], parent=None, counts=counts,
                               lda=0.0, cls=None, color='#BBBBBB', xy=(0.0, 0.0))

    # ---------------------------------------------------------------
    # Parent-child inference via count correlation
    # For each node at rank[i], find the parent at rank[i-1] whose
    # count vector has the highest Pearson r with the child's vector.
    # Fall back to first node of parent rank if no good correlation.
    # ---------------------------------------------------------------
    for r_idx in range(len(present_ranks)):
        rank     = present_ranks[r_idx]
        children = [nid for nid, n in nodes.items() if n['rank'] == rank]

        if r_idx == 0:
            parent_pool = ['__root__']
        else:
            parent_rank = present_ranks[r_idx - 1]
            parent_pool = [nid for nid, n in nodes.items()
                           if n['rank'] == parent_rank]

        if not parent_pool:
            parent_pool = ['__root__']

        # Pre-stack parent counts
        par_counts = np.vstack([nodes[p]['counts'] for p in parent_pool])  # (P, S)

        for child_id in children:
            cv = nodes[child_id]['counts']
            if cv.sum() == 0 or par_counts.shape[1] == 0:
                best_par = parent_pool[0]
            else:
                # Pearson r between child and each parent
                corrs = []
                for pv in par_counts:
                    if pv.std() < 1e-12 or cv.std() < 1e-12:
                        corrs.append(0.0)
                    else:
                        corrs.append(float(np.corrcoef(cv, pv)[0, 1]))
                best_par = parent_pool[int(np.argmax(corrs))]

            nodes[child_id]['parent'] = best_par
            nodes[best_par]['children'].append(child_id)

    return nodes, present_ranks


def _assign_angles(node_id, nodes, start, end):
    """Recursively assign angular sectors to all nodes (DFS)."""
    node     = nodes[node_id]
    node['angle_start'] = start
    node['angle_end']   = end
    node['angle_mid']   = (start + end) / 2.0
    children = node['children']
    if not children:
        return
    span = end - start
    for ch_id in children:
        nodes[ch_id]['_w'] = max(1, len(_leaf_count(ch_id, nodes)))
    total_w = sum(nodes[c]['_w'] for c in children)
    cur = start
    for ch_id in children:
        w     = nodes[ch_id]['_w']
        share = span * w / total_w
        _assign_angles(ch_id, nodes, cur, cur + share)
        cur += share


def _leaf_count(node_id, nodes):
    """Return list of leaf descendants (nodes with no children)."""
    ch = nodes[node_id]['children']
    if not ch:
        return [node_id]
    leaves = []
    for c in ch:
        leaves.extend(_leaf_count(c, nodes))
    return leaves


def generate_cladogram(nodes, present_ranks, biomarkers_df, color_map,
                        output_base, fmt, title_suffix=''):
    """
    Draw the classic LEfSe circular cladogram.

    Parameters
    ----------
    nodes          : tree built by build_taxonomy_tree()
    present_ranks  : ordered list of rank strings
    biomarkers_df  : DataFrame with columns Feature, Class, Abs_LDA
    color_map      : {class_name: hex_color}
    output_base    : path without extension
    fmt            : 'png'|'pdf'|'tiff'
    """
    print(f"    [*] Drawing cladogram → {output_base}_cladogram")

    # Annotate significant nodes
    sig_lookup = {}   # name_lower → (cls, abs_lda)
    for _, row in biomarkers_df.iterrows():
        sig_lookup[row['Feature'].strip().lower()] = (row['Class'], row['Abs_LDA'])

    n_marked = 0
    for nid, node in nodes.items():
        key = node['name'].strip().lower()
        if key in sig_lookup:
            cls, lda = sig_lookup[key]
            node['cls']   = cls
            node['color'] = color_map.get(cls, '#888888')
            node['lda']   = lda
            n_marked += 1
    print(f"    [*] Nodes marked as significant: {n_marked}")

    # Assign angles to entire tree
    _assign_angles('__root__', nodes, 0, 2 * np.pi)

    # Radii per rank level
    N_RINGS      = len(present_ranks)
    R_STEP       = 1.0 / (N_RINGS + 1)       # gap between rings
    rank_radius  = {r: (i + 1) * R_STEP for i, r in enumerate(present_ranks)}
    ROOT_R       = 0.0

    # ----------------------------------------------------------------
    # Figure
    # ----------------------------------------------------------------
    fig_size = max(14, N_RINGS * 2.5 + 4)
    fig, ax  = plt.subplots(figsize=(fig_size, fig_size))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.set_aspect('equal')
    ax.axis('off')

    # Draw faint ring guide circles
    for r in rank_radius.values():
        circle = plt.Circle((0, 0), r, color='#EEEEEE',
                             fill=False, linewidth=0.6, zorder=1)
        ax.add_patch(circle)

    # ----------------------------------------------------------------
    # Draw edges and nodes
    # ----------------------------------------------------------------
    MAX_LDA      = max((n.get('lda', 0) for n in nodes.values() if n.get('lda', 0) > 0), default=1.0)
    NODE_MIN     = 30    # pt² for non-significant
    NODE_MAX     = 350   # pt² for max-LDA significant

    drawn_labels = []  # (x, y, text, color) — for collision-aware labels

    for nid, node in nodes.items():
        if nid == '__root__':
            continue
        r   = rank_radius[node['rank']]
        ang = node['angle_mid']
        x   = r * np.cos(ang)
        y   = r * np.sin(ang)
        node['xy'] = (x, y)

        # Edge to parent
        par_id = node['parent']
        if par_id and par_id in nodes:
            par = nodes[par_id]
            if par_id == '__root__':
                px, py = 0.0, 0.0
            else:
                px, py = par.get('xy', (0.0, 0.0))
            ax.plot([px, x], [py, y],
                    color='#CCCCCC', linewidth=0.7, zorder=2, solid_capstyle='round')

        # Node circle
        is_sig  = node['cls'] is not None
        color   = node['color']
        size    = NODE_MIN
        if is_sig and MAX_LDA > 0:
            size = NODE_MIN + (NODE_MAX - NODE_MIN) * (node['lda'] / MAX_LDA)
        edgec   = '#555555' if is_sig else '#AAAAAA'
        lw      = 1.2      if is_sig else 0.5
        ax.scatter(x, y, s=size, c=color, edgecolors=edgec,
                   linewidths=lw, zorder=4)

        # Collect labels for significant nodes
        if is_sig:
            drawn_labels.append((x, y, node['name'], color, node['lda']))

    # ----------------------------------------------------------------
    # Labels for significant nodes — radial placement with collision avoidance
    # ----------------------------------------------------------------
    LABEL_R_PAD = 0.06   # extra radial push beyond node centre
    FONT_SIZE   = 7.5

    # Sort by angle to make collision detection meaningful
    drawn_labels.sort(key=lambda t: np.arctan2(t[1], t[0]))

    placed = []  # (x_text, y_text)

    for (nx, ny, name, color, lda) in drawn_labels:
        dist  = np.hypot(nx, ny)
        ang   = np.arctan2(ny, nx)
        r_txt = dist + LABEL_R_PAD

        tx = r_txt * np.cos(ang)
        ty = r_txt * np.sin(ang)

        # Push away from already-placed labels
        for attempt in range(20):
            too_close = any(np.hypot(tx - px, ty - py) < 0.07 for px, py in placed)
            if not too_close:
                break
            r_txt += 0.05
            tx = r_txt * np.cos(ang)
            ty = r_txt * np.sin(ang)

        placed.append((tx, ty))

        ha = 'left' if tx >= 0 else 'right'
        rotation = np.degrees(ang)
        if tx < 0:
            rotation += 180

        ax.text(tx, ty, name,
                ha=ha, va='center',
                fontsize=FONT_SIZE, fontweight='bold',
                color=color, rotation=rotation,
                rotation_mode='anchor', zorder=5)

        # Leader dot-line from node to text
        ax.plot([nx, tx], [ny, ty],
                color=color, linewidth=0.5, linestyle=':', zorder=3, alpha=0.6)

    # ----------------------------------------------------------------
    # Ring level labels (rank names on the right side)
    # ----------------------------------------------------------------
    for rank, r in rank_radius.items():
        ax.text(r * 1.01, -0.015, rank.capitalize(),
                fontsize=7, color='#999999', ha='left', va='top',
                style='italic', zorder=6)

    # ----------------------------------------------------------------
    # Legend
    # ----------------------------------------------------------------
    patches = [mpatches.Patch(color=c, label=g) for g, c in color_map.items()]
    ax.legend(handles=patches, title='Class', title_fontsize=9,
              fontsize=8.5, loc='lower left',
              bbox_to_anchor=(0.01, 0.01),
              frameon=True, framealpha=0.9, edgecolor='#CCCCCC')

    # Size legend (LDA scale)
    for lda_val, label in [(2, '|LDA|=2'), (3, '|LDA|=3'), (4, '|LDA|=4')]:
        if lda_val <= MAX_LDA + 0.5:
            size = NODE_MIN + (NODE_MAX - NODE_MIN) * (lda_val / MAX_LDA)
            ax.scatter([], [], s=size, c='#888888', edgecolors='#555555',
                       linewidths=1.0, label=label, zorder=0)
    ax.legend(handles=patches + [
        mpatches.Patch(color='w', label='─── Node size = |LDA|')
    ], title='Class', title_fontsize=9, fontsize=8.5,
              loc='lower left', bbox_to_anchor=(0.01, 0.01),
              frameon=True, framealpha=0.9, edgecolor='#CCCCCC')

    title = f'LEfSe Cladogram{title_suffix}'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=16)

    lim = 1.35
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    plt.tight_layout(pad=1.0)
    cladogram_base = f"{output_base}_cladogram"
    export_topological_projection(fig, cladogram_base, fmt)
    print(f"    [*] Cladogram saved: {cladogram_base}.{fmt}")
    plt.close(fig)


# ===========================================================================
# Main LEfSe pipeline
# ===========================================================================

def run_lefse(df_rank, rank_level, metadata_path, category_col, sample_id_col,
              output_base, fmt,
              lda_threshold=2.0, kw_alpha=0.05, wilcox_alpha=0.05,
              top_n=30, sort_by_lda=False, label_col=None, no_table=False,
              df_full=None, no_cladogram=False):

    print(f"[*] Running LEfSe → {output_base}")

    # -------------------------------------------------------------------
    # 1. Resolve feature labels
    # -------------------------------------------------------------------
    df_work = df_rank.copy()
    if label_col and label_col in df_work.columns:
        print(f"    [*] Feature labels from column: '{label_col}'")
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

    # -------------------------------------------------------------------
    # 2. Load metadata & align samples
    # -------------------------------------------------------------------
    if metadata_path.endswith('.csv'):
        meta = pd.read_csv(metadata_path)
    else:
        meta = pd.read_excel(metadata_path)

    # ------------------------------------------------------------------
    # Sample-name normalisation — mirrors 05_generate_PCoA_PieChart.py
    # Strips lane/run suffixes (e.g. _L6, _L8) by splitting on '_' and
    # comparing only the first token, case-insensitively.
    # ------------------------------------------------------------------
    meta[sample_id_col] = meta[sample_id_col].astype(str).str.strip().str.lower()
    meta_map = dict(zip(meta[sample_id_col], meta[category_col].astype(str)))

    all_samples   = list(df_counts.columns)
    meta_keys     = list(meta_map.keys())   # already lowercased
    groups_raw    = []
    for name in all_samples:
        clean = name.split('_')[0].strip().lower()
        val   = meta_map.get(clean, 'Unknown')

        # Fallback: longest metadata key that is a prefix of `clean`
        if val in ('Unknown', 'nan'):
            candidates = [k for k in meta_keys if clean.startswith(k)]
            if candidates:
                best_key = max(candidates, key=len)
                val      = meta_map[best_key]
                print(f"    [~] Prefix match: '{name}' → '{best_key}' (group: {val})")

        if val == 'nan':
            val = 'Unknown'
        if val == 'Unknown':
            print(f"    [!] Topological Mismatch: '{name}' (cleaned: '{clean}') NOT FOUND in metadata.")
        groups_raw.append(val)

    # Drop unmatched samples
    keep_mask      = [g != 'Unknown' for g in groups_raw]
    dropped        = [s for s, k in zip(all_samples, keep_mask) if not k]
    common_samples = [s for s, k in zip(all_samples, keep_mask) if k]
    groups         = [g for g, k in zip(groups_raw,  keep_mask) if k]

    if dropped:
        print(f"    [!] Dropped {len(dropped)} unmatched sample(s): {dropped}")
    if len(common_samples) < 3:
        raise ValueError("CRITICAL: fewer than 3 samples matched metadata — cannot run LEfSe.")

    df_counts = df_counts[common_samples]

    unique_groups = sorted(set(groups))
    n_groups      = len(unique_groups)
    print(f"    [*] Groups ({n_groups}): {unique_groups}")
    print(f"    [*] Samples matched: {len(common_samples)}")

    # -------------------------------------------------------------------
    # 3. Relative abundance & filter zero-variance features
    # -------------------------------------------------------------------
    ra_df = relative_abundance(df_counts)
    ra_df = ra_df.fillna(0)
    ra_df = ra_df.loc[ra_df.var(axis=1) > 0]
    print(f"    [*] Features after zero-variance filter: {len(ra_df)}")

    # -------------------------------------------------------------------
    # 4. Kruskal-Wallis
    # -------------------------------------------------------------------
    kw_results = kruskal_wallis_test(ra_df, groups, alpha=kw_alpha)
    print(f"    [*] Features passing Kruskal-Wallis (p<{kw_alpha}): {len(kw_results)}")
    if not kw_results:
        print("    [!] No features passed Kruskal-Wallis. Exiting.")
        return

    # -------------------------------------------------------------------
    # 5. Pairwise Wilcoxon
    # -------------------------------------------------------------------
    winning_classes = wilcoxon_pairwise_test(
        ra_df, groups, list(kw_results.keys()), alpha=wilcox_alpha
    )
    print(f"    [*] Features passing Wilcoxon (p<{wilcox_alpha}): {len(winning_classes)}")
    if not winning_classes:
        print("    [!] No features passed Wilcoxon. Exiting.")
        return

    # -------------------------------------------------------------------
    # 6. LDA scores
    # -------------------------------------------------------------------
    lda_scores = compute_lda_scores(ra_df, groups, list(winning_classes.keys()), winning_classes)

    # -------------------------------------------------------------------
    # 7. Build results table & apply LDA threshold
    # -------------------------------------------------------------------
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
            'LDA_score':     round(lda,    4),
            'Abs_LDA':       round(abs_lda,4),
            'KW_H':          round(H,      4),
            'KW_p':          round(kw_p,   6),
            'Mean_RA':       round(float(ra_df.loc[feat].mean()), 6),
        })

    if not records:
        print(f"    [!] No features passed |LDA| >= {lda_threshold}. Try lowering --lda_threshold.")
        return

    df_results = pd.DataFrame(records).sort_values('Abs_LDA', ascending=False)
    print(f"    [*] Significant biomarkers (|LDA| >= {lda_threshold}): {len(df_results)}")

    # Subset to top_n
    df_plot = df_results.head(top_n).copy()

    # -------------------------------------------------------------------
    # 8. Group-level mean RA table (extra sheet)
    # -------------------------------------------------------------------
    group_means = []
    for feat in df_plot['Feature']:
        row = {'Feature': feat}
        for ug in unique_groups:
            samp = [s for s, g in zip(common_samples, groups) if g == ug]
            row[f'MeanRA_{ug}'] = round(float(ra_df.loc[feat, samp].mean()), 6)
        group_means.append(row)
    df_group_means = pd.DataFrame(group_means)

    # -------------------------------------------------------------------
    # 9. Plot — horizontal bar chart
    # -------------------------------------------------------------------
    color_map = {g: COLORS[i % len(COLORS)] for i, g in enumerate(unique_groups)}

    if sort_by_lda:
        df_plot = df_plot.sort_values('LDA_score')
    else:
        # Group bars by winning class, then by LDA within class
        df_plot = df_plot.sort_values(['Class', 'Abs_LDA'], ascending=[True, False])

    n_bars   = len(df_plot)
    fig_h    = max(8, n_bars * 0.38 + 2)
    fig, ax  = plt.subplots(figsize=(13, fig_h))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('#F7F7F7')

    # Draw bars
    for i, row in enumerate(df_plot.itertuples()):
        lda_val = row.LDA_score
        color   = color_map[row.Class]
        bar_val = lda_val if sort_by_lda else row.Abs_LDA
        bar_dir = 1 if lda_val >= 0 else -1

        if not sort_by_lda:
            ax.barh(i, bar_val, color=color, edgecolor='white', linewidth=0.5,
                    height=0.72, zorder=3)
        else:
            ax.barh(i, lda_val, color=color, edgecolor='white', linewidth=0.5,
                    height=0.72, zorder=3)

        # Score label inside/outside bar
        label_x = bar_val + 0.05 if not sort_by_lda else (lda_val + 0.05 * bar_dir)
        ax.text(label_x, i, f"{row.Abs_LDA:.2f}",
                va='center', ha='left', fontsize=8.5, color='#333333', zorder=4)

    # Feature labels on y-axis
    ax.set_yticks(range(n_bars))
    ax.set_yticklabels(df_plot['Feature'].tolist(), fontsize=9.5)

    # Axis styling
    ax.axvline(0, color='#888888', linewidth=0.8, linestyle='--', zorder=2)
    ax.set_xlabel('LDA Score (log₁₀)', fontsize=12, labelpad=8)
    ax.set_title(
        f'LEfSe — {rank_level.capitalize()} level   '
        f'[|LDA| ≥ {lda_threshold} | KW p<{kw_alpha} | Wilcoxon p<{wilcox_alpha}]',
        fontsize=13, fontweight='bold', pad=14
    )
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#CCCCCC')
    ax.tick_params(axis='both', which='both', length=0)
    ax.grid(axis='x', color='white', linewidth=1.2, zorder=1)
    ax.set_xlim(left=0 if not sort_by_lda else None)

    # Legend — one patch per group
    patches = [mpatches.Patch(color=color_map[g], label=g) for g in unique_groups]
    ax.legend(handles=patches, title='Class', title_fontsize=10,
              fontsize=9.5, loc='lower right',
              frameon=True, framealpha=0.9, edgecolor='#CCCCCC')

    # Subtitle: n biomarkers found
    ax.text(0.01, 1.01, f"n = {len(df_results)} biomarkers | showing top {n_bars}",
            transform=ax.transAxes, fontsize=9, color='#666666', ha='left')

    plt.tight_layout(pad=1.5)

    export_topological_projection(fig, output_base, fmt)
    print(f"    [*] Plot saved: {output_base}.{fmt}")
    plt.close(fig)

    # -------------------------------------------------------------------
    # 10. Export Excel
    # -------------------------------------------------------------------
    if not no_table:
        export_summary_table(
            df_results, output_base,
            sheet_name='LEfSe_Biomarkers',
            extra_sheets={
                'Group_Mean_RelAbund': df_group_means,
            }
        )

    # -------------------------------------------------------------------
    # 11. Cladogram
    # -------------------------------------------------------------------
    if not no_cladogram:
        src = df_full if df_full is not None else df_rank
        nodes, present_ranks = build_taxonomy_tree(src, common_samples, label_col=label_col)
        if present_ranks:
            generate_cladogram(
                nodes, present_ranks, df_results, color_map,
                output_base, fmt,
                title_suffix=f' — {rank_level.capitalize()} level'
            )
        else:
            print("    [!] Could not build taxonomy tree — skipping cladogram.")


# ===========================================================================
# Entry point
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description='LEfSe — Linear Discriminant Analysis Effect Size for Kraken2 pipeline output.'
    )
    parser.add_argument('-d',  '--data',          type=str, required=True,
                        help='Input Excel file (taxonomic_classification_clean.xlsx).')
    parser.add_argument('-r',  '--rank',           type=str, required=True,
                        help='Taxonomic rank to analyse (e.g. genus, species, domain).')
    parser.add_argument('-m',  '--metadata',       type=str, required=True,
                        help='Metadata CSV or Excel with sample IDs and group column.')
    parser.add_argument('-c',  '--category',       type=str, required=True,
                        help='Metadata column to group samples by.')
    parser.add_argument('-id', '--sample_id',      type=str, default='SampleID',
                        help='Metadata column for sample IDs (default: SampleID).')
    parser.add_argument('-fmt','--format',         type=str,
                        choices=['pdf', 'png', 'tiff'], default='pdf',
                        help='Output format: pdf, png, or tiff (default: pdf).')
    parser.add_argument('-o',  '--output',         type=str, default=None,
                        help='Custom output base name (no extension).')
    parser.add_argument('-org','--organism',       type=str, default='Microbiome',
                        help='Label prefix for output files.')
    parser.add_argument('--lda_threshold',         type=float, default=2.0,
                        help='Minimum |LDA| score cut-off (default: 2.0).')
    parser.add_argument('--kw_alpha',              type=float, default=0.05,
                        help='Kruskal-Wallis significance cut-off (default: 0.05).')
    parser.add_argument('--wilcox_alpha',          type=float, default=0.05,
                        help='Wilcoxon rank-sum significance cut-off (default: 0.05).')
    parser.add_argument('--top',                   type=int,   default=30,
                        help='Maximum number of biomarkers to display (default: 30).')
    parser.add_argument('--rank_by_lda',           action='store_true',
                        help='Sort bars by signed LDA score instead of grouping by class.')
    parser.add_argument('--label_col',             type=str, default=None,
                        help='Column to use as feature label (e.g. "Name", "Scientific Name").')
    parser.add_argument('--no_cladogram',          action='store_true',
                        help='Skip cladogram plot.')
    parser.add_argument('--no_table',              action='store_true',
                        help='Skip Excel table export.')

    args = parser.parse_args()

    OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'results', 'LEfSe')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    safe_org  = args.organism.replace(' ', '_')
    base_name = args.output if args.output else f"{safe_org}_{args.rank}_LEfSe"
    output_base = os.path.join(OUTPUT_DIR, base_name)

    df      = pd.read_excel(args.data, sheet_name=0)
    df_rank = df[df['Rank'].str.lower() == args.rank.lower()].copy()

    if df_rank.empty:
        raise ValueError(f"No data found for rank: {args.rank}")

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
        df_full       = df,
        no_cladogram  = args.no_cladogram,
    )


if __name__ == '__main__':
    main()