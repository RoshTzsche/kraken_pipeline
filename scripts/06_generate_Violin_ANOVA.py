import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from itertools import combinations
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
import math

# ---------------------------------------------------------------------------
# Unified color palette — identical across all pipeline scripts (04, 05, 07)
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
# Global typography — mirrors R Tahoma/bold preference
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

# Strip header constants (mirrors ggplot2 facet strip defaults)
STRIP_BG   = "#EBEBEB"
STRIP_EDGE = "#999999"
STRIP_H_IN = 0.30          # strip height in inches


# ---------------------------------------------------------------------------
# ggplot2-style discrete hue palette (kept for backward-compatibility with taxa mode)
# ---------------------------------------------------------------------------
def ggplot2_palette(n):
    import colorsys
    return [colorsys.hls_to_rgb(h / n, 0.65, 0.70) for h in range(n)]


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
def export_figure(fig, output_base, fmt, pad_inches=0.35):
    """
    Projects the continuous mathematical representation (matplotlib Figure) into
    a defined discrete or continuous state space (PNG, TIFF, or PDF).
    Maintains 300 DPI high-frequency spatial resolution for clinical rasters.
    """
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
    """
    Serializes the computed summary DataFrame into a formatted Excel workbook
    with auto-adjusted column widths for readability.
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


# ---------------------------------------------------------------------------
# Compact Letter Display — Tukey HSD (parametric, for taxa mode)
# ---------------------------------------------------------------------------
def calculate_cld(df, val_col, group_col):
    mean_series = df.groupby(group_col)[val_col].mean().sort_values(ascending=False)
    groups      = mean_series.index.tolist()

    if len(groups) < 2:
        return {g: "a" for g in groups}

    group_data = [df[df[group_col] == g][val_col].dropna().values for g in groups]
    group_data = [g for g in group_data if len(g) > 0]

    if len(group_data) < 2 or df[val_col].nunique() <= 1:
        return {g: "a" for g in groups}

    _, p_val = stats.f_oneway(*group_data)
    if p_val > 0.05 or pd.isna(p_val):
        return {g: "a" for g in groups}

    tukey      = pairwise_tukeyhsd(endog=df[val_col], groups=df[group_col], alpha=0.05)
    sig_matrix = pd.DataFrame(False, index=groups, columns=groups)
    for row in tukey.summary().data[1:]:
        g1, g2, reject = row[0], row[1], row[-1]
        if str(reject).strip().lower() == "true":
            sig_matrix.loc[g1, g2] = True
            sig_matrix.loc[g2, g1] = True

    letters      = {g: "" for g in groups}
    current_char = ord("a")
    for g1 in groups:
        if not letters[g1]:
            letter = chr(current_char)
            current_char += 1
            letters[g1] += letter
            for g2 in groups:
                if g1 != g2 and not sig_matrix.loc[g1, g2]:
                    if letter not in letters[g2]:
                        letters[g2] += letter

    return letters if any(letters.values()) else {g: "a" for g in groups}


# ---------------------------------------------------------------------------
# Compact Letter Display — Kruskal-Wallis + Mann-Whitney Bonferroni
# (non-parametric, for alpha diversity mode)
# ---------------------------------------------------------------------------
def compute_kruskal_cld(df, val_col, group_col):
    """
    Non-parametric CLD via Kruskal-Wallis omnibus test followed by pairwise
    Mann-Whitney U comparisons with Bonferroni correction.
    Groups sharing a letter are NOT significantly different from each other.
    """
    # Sort groups by mean value descending (matches CLD convention)
    mean_vals = df.groupby(group_col)[val_col].mean().sort_values(ascending=False)
    groups    = mean_vals.index.tolist()

    if len(groups) < 2:
        return {g: 'a' for g in groups}

    group_data  = [df[df[group_col] == g][val_col].dropna().values for g in groups]
    valid_data  = [d for d in group_data if len(d) > 0]

    if len(valid_data) < 2 or df[val_col].nunique() <= 1:
        return {g: 'a' for g in groups}

    # Kruskal-Wallis omnibus test
    try:
        _, kw_p = stats.kruskal(*valid_data)
    except Exception:
        return {g: 'a' for g in groups}

    if kw_p > 0.05 or np.isnan(kw_p):
        return {g: 'a' for g in groups}

    # Pairwise Mann-Whitney U with Bonferroni correction
    pair_list = list(combinations(range(len(groups)), 2))
    p_values  = []
    for i, j in pair_list:
        d1, d2 = group_data[i], group_data[j]
        if len(d1) < 1 or len(d2) < 1:
            p_values.append(1.0)
        else:
            try:
                _, p = stats.mannwhitneyu(d1, d2, alternative='two-sided')
            except Exception:
                p = 1.0
            p_values.append(p)

    if not p_values:
        return {g: 'a' for g in groups}

    reject, _, _, _ = multipletests(p_values, alpha=0.05, method='bonferroni')

    # Build significance matrix
    sig_mat = {g: {g2: False for g2 in groups} for g in groups}
    for (i, j), rej in zip(pair_list, reject):
        if rej:
            sig_mat[groups[i]][groups[j]] = True
            sig_mat[groups[j]][groups[i]] = True

    # Assign CLD letters (same algorithm as calculate_cld for visual consistency)
    letters      = {g: '' for g in groups}
    current_char = ord('a')
    for g1 in groups:
        if not letters[g1]:
            letter = chr(current_char)
            current_char += 1
            letters[g1] += letter
            for g2 in groups:
                if g1 != g2 and not sig_mat[g1][g2]:
                    if letter not in letters[g2]:
                        letters[g2] += letter

    return letters if any(letters.values()) else {g: 'a' for g in groups}


# ---------------------------------------------------------------------------
# Alpha Diversity Index Computation
# ---------------------------------------------------------------------------
def compute_alpha_diversity(counts_series):
    """
    Computes Shannon H', Simpson (1-D), and bias-corrected Chao1 from a
    pandas Series of raw integer read counts per taxon for a single sample.
    """
    counts = counts_series.values.astype(float)
    counts = counts[counts > 0]       # exclude zero-count taxa
    N      = counts.sum()

    if N == 0:
        return {'Shannon': 0.0, 'Simpson': 0.0, 'Chao1': 0.0}

    # Shannon H'
    p       = counts / N
    shannon = float(-np.sum(p * np.log(p + 1e-15)))

    # Simpson (1-D) — bias-corrected form
    if N > 1:
        D       = np.sum(counts * (counts - 1)) / (N * (N - 1))
        simpson = float(1.0 - D)
    else:
        simpson = 0.0

    # Chao1 bias-corrected
    s_obs = len(counts)
    f1    = int(np.sum(counts == 1))   # singletons
    f2    = int(np.sum(counts == 2))   # doubletons
    chao1 = float(s_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1)))

    return {'Shannon': shannon, 'Simpson': simpson, 'Chao1': chao1}


# ---------------------------------------------------------------------------
# theme_bw — all 4 spines visible ("caged" / enclosed panels)
# ---------------------------------------------------------------------------
def apply_theme_bw(ax):
    """
    Mirrors ggplot2 theme_bw():
      - White panel background, no grid
      - All four spines visible and black (full enclosure)
      - Ticks outward
    """
    ax.set_facecolor("white")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")
    ax.tick_params(axis="both", which="both", direction="out",
                   colors="black", width=0.8, length=4)
    ax.tick_params(axis="x", labelcolor="black", labelsize=9, rotation=45)
    ax.tick_params(axis="y", labelcolor="black", labelsize=9)


# ---------------------------------------------------------------------------
# Draw one violin facet panel (used by taxa mode)
# No per-panel title or axis labels — those are shared at figure level
# ---------------------------------------------------------------------------
def draw_violin_panel(ax, sub_df, group_order, palette_map):
    groups        = group_order
    data_by_group = [sub_df[sub_df["Group"] == g]["Abundance"].dropna().values
                     for g in groups]

    # Violin bodies — only for groups with >1 data point
    valid_pos  = [i for i, d in enumerate(data_by_group) if len(d) > 1]
    valid_data = [data_by_group[i] for i in valid_pos]

    if valid_data:
        parts = ax.violinplot(valid_data, positions=valid_pos,
                              widths=0.70, showmeans=False,
                              showmedians=False, showextrema=False)
        for idx, body in zip(valid_pos, parts["bodies"]):
            body.set_facecolor(palette_map[groups[idx]])
            body.set_edgecolor("black")
            body.set_linewidth(1.1)
            body.set_alpha(0.60)

    # Mean point + SE error bar
    for i, g in enumerate(groups):
        vals = data_by_group[i]
        if len(vals) == 0:
            continue
        mean_val = np.mean(vals)
        se_val   = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
        ax.errorbar(i, mean_val, yerr=se_val,
                    fmt="none", ecolor="black",
                    elinewidth=1.1, capsize=4, capthick=1.1, zorder=4)
        ax.plot(i, mean_val, "o", color="black",
                markersize=4.5, zorder=5,
                markeredgecolor="black", markeredgewidth=0.7)

    # CLD letters
    try:
        letters_dict = calculate_cld(sub_df, "Abundance", "Group")
    except Exception:
        letters_dict = {g: "" for g in groups}

    y_vals = sub_df["Abundance"].dropna()
    y_max  = float(y_vals.max()) if len(y_vals) else 1.0
    y_pos  = y_max * 1.12
    ax.set_ylim(bottom=0, top=y_pos * 1.12)

    for i, g in enumerate(groups):
        ax.text(i, y_pos, letters_dict.get(g, ""),
                ha="center", va="bottom",
                fontsize=11, fontweight="bold", color="black")

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, ha="right")
    apply_theme_bw(ax)


# ---------------------------------------------------------------------------
# Draw facet strip headers AFTER tight_layout
# ---------------------------------------------------------------------------
def draw_strip_headers(fig, axes_flat, taxa_list):
    """
    Draws a full-width gray rectangle + italic bold taxon name above each
    panel, mimicking ggplot2 facet_wrap strip headers.
    Must be called after tight_layout so ax.get_position() is finalised.
    """
    fig_h_in     = fig.get_figheight()
    strip_h_frac = STRIP_H_IN / fig_h_in

    for i, taxon in enumerate(taxa_list):
        ax  = axes_flat[i]
        pos = ax.get_position()      # Bbox in figure-fraction coords

        # Gray background strip
        rect = Rectangle(
            (pos.x0, pos.y1),
            pos.width, strip_h_frac,
            transform=fig.transFigure,
            facecolor=STRIP_BG,
            edgecolor=STRIP_EDGE,
            linewidth=0.8,
            clip_on=False,
            zorder=5,
        )
        fig.add_artist(rect)

        # Taxon label centred on the strip
        fig.text(
            pos.x0 + pos.width / 2,
            pos.y1 + strip_h_frac / 2,
            taxon,
            ha="center", va="center",
            fontsize=10, fontweight="bold", fontstyle="italic",
            color="black", clip_on=False, zorder=6,
        )

        # Align top spine colour with strip edge
        ax.spines["top"].set_color(STRIP_EDGE)


# ---------------------------------------------------------------------------
# Main pipeline — TAXA MODE (original behaviour, ANOVA/Tukey CLD)
# ---------------------------------------------------------------------------
def generate_taxa_violin_plots(data_path, metadata_path, category_col,
                                sample_id_col, rank_level, threshold,
                                organism_name, output_base, fmt, no_table=False):
    """
    Generates per-taxon violin plots (ANOVA + Tukey HSD CLD) in ggplot2 facet_wrap style.
    Optionally exports a long-format abundance summary table (.xlsx).
    """
    print(f"[*] Generating statistical violin plots (ANOVA/Tukey) → {output_base}")

    # 1. Load & filter
    df      = pd.read_excel(data_path, sheet_name=0)
    df_rank = df[df["Rank"].str.lower() == rank_level.lower()].copy()
    if df_rank.empty:
        raise ValueError(f"No data found for taxonomic level: {rank_level}")

    meta_cols   = {"Rank", "TaxID", "original_header", "Name", "Scientific Name"}
    sample_cols = [c for c in df_rank.columns if c not in meta_cols]
    df_counts   = df_rank[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    # 2. Relative abundance
    rel_abund          = df_counts.div(df_counts.sum(axis=0), axis=1) * 100
    rel_abund["Taxon"] = df_rank["Name"].values

    # 3. Threshold filter
    mean_rel    = rel_abund[sample_cols].mean(axis=1)
    top_taxa_df = rel_abund[mean_rel >= threshold * 100]
    if top_taxa_df.empty:
        print(f"[!] No taxa passed the {threshold*100:.2f}% threshold.")
        return

    # 4. Melt
    melted = top_taxa_df.melt(id_vars=["Taxon"], value_vars=sample_cols,
                               var_name="Sample", value_name="Abundance")
    melted["CleanSample"] = melted["Sample"].apply(
        lambda x: str(x).split("_")[0].strip().lower())

    # 5. Metadata
    meta_df = (pd.read_csv(metadata_path) if metadata_path.lower().endswith(".csv")
               else pd.read_excel(metadata_path))
    meta_df[sample_id_col] = meta_df[sample_id_col].astype(str).str.strip().str.lower()

    if pd.api.types.is_numeric_dtype(meta_df[category_col]):
        print(f"    [*] Continuous variable '{category_col}' — Quantile Binning.")
        meta_df[category_col] = pd.qcut(meta_df[category_col].dropna(), q=4).astype(str)

    meta_map        = dict(zip(meta_df[sample_id_col],
                               meta_df[category_col].astype(str)))
    melted["Group"] = melted["CleanSample"].map(meta_map)
    melted          = melted.dropna(subset=["Group"])
    if melted.empty:
        raise ValueError("No Sample IDs matched between abundance data and metadata.")

    # 6. Group order + palette
    group_order = sorted(melted["Group"].unique())
    palette_map = dict(zip(group_order, ggplot2_palette(len(group_order))))

    # 7. Figure dimensions
    taxa_list  = list(melted["Taxon"].unique())
    n_vars     = len(taxa_list)
    ncols      = min(3, n_vars)
    nrows      = math.ceil(n_vars / ncols)

    legend_w_in = 1.8
    fig_w       = ncols * 4.8 + legend_w_in
    fig_h       = nrows * 5.2 + 0.7

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), facecolor="white")
    fig.patch.set_facecolor("white")

    if n_vars == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # 8. Draw each panel
    for i, taxon in enumerate(taxa_list):
        sub_df = melted[melted["Taxon"] == taxon].dropna(subset=["Abundance"])
        draw_violin_panel(axes[i], sub_df, group_order, palette_map)

    # Hide unused axes
    for j in range(n_vars, len(axes)):
        axes[j].set_visible(False)

    # X tick labels only on bottom-most panel of each column
    for col in range(ncols):
        col_panel_indices = [
            col + row * ncols
            for row in range(nrows)
            if col + row * ncols < n_vars
        ]
        for idx in col_panel_indices[:-1]:
            axes[idx].tick_params(axis="x", labelbottom=False)
            axes[idx].set_xticklabels([])
        axes[col_panel_indices[-1]].tick_params(
            axis="x", labelbottom=True, labelrotation=45, labelsize=9
        )

    # 9. tight_layout with reserved margins
    legend_frac = legend_w_in / fig_w
    strip_frac  = (STRIP_H_IN * nrows + 0.3) / fig_h

    fig.tight_layout(
        rect=[0.04, 0.06, 1.0 - legend_frac - 0.01, 1.0 - strip_frac],
        pad=1.6, h_pad=4.0, w_pad=2.0,
    )

    # 10. Shared axis labels
    fig.supxlabel(category_col, fontsize=12, fontweight="bold", y=0.01, va="bottom")
    fig.supylabel(f"{organism_name} Rel. Abundance (%)",
                  fontsize=11, fontweight="bold", x=0.01, va="center")

    # 11. Right-side vertical legend
    legend_patches = [
        mpatches.Patch(facecolor=palette_map[g], edgecolor="black",
                       linewidth=0.7, alpha=0.7, label=g)
        for g in group_order
    ]
    fig.legend(
        handles=legend_patches,
        title=category_col,
        title_fontsize=9,
        fontsize=9,
        loc="center left",
        bbox_to_anchor=(1.0 - legend_frac + 0.01, 0.5),
        frameon=True,
        framealpha=0.95,
        edgecolor="#CCCCCC",
        ncol=1,
    )

    # 12. Facet strip headers (drawn AFTER tight_layout — positions are final)
    draw_strip_headers(fig, axes, taxa_list)

    # 13. Export
    export_figure(fig, output_base, fmt)
    plt.close(fig)
    print(f"    [✓] Plot saved: {output_base}.{fmt}")

    # 14. Long-format summary table
    if not no_table:
        table_df = (
            melted[['Taxon', 'Sample', 'Group', 'Abundance']]
            .rename(columns={'Abundance': 'Rel_Abundance_pct'})
            .round({'Rel_Abundance_pct': 4})
            .sort_values(['Taxon', 'Group', 'Sample'])
            .reset_index(drop=True)
        )
        export_summary_table(table_df, output_base, sheet_name='Taxa_Abundance')


# ---------------------------------------------------------------------------
# Main pipeline — ALPHA DIVERSITY MODE (Shannon, Simpson, Chao1)
# ---------------------------------------------------------------------------
def generate_alpha_diversity_violins(data_path, metadata_path, category_col,
                                      sample_id_col, rank_level, organism_name,
                                      output_base, fmt, no_table=False):
    """
    Generates a 1×3 multi-panel violin figure for three alpha diversity indices:
      - Shannon
      - Simpson
      - Chao1

    Violins are grouped by a metadata category column. CLD letters are computed
    via Kruskal-Wallis + pairwise Mann-Whitney U (Bonferroni correction) and
    placed above each violin. 
    Optionally exports a per-sample summary table (.xlsx).
    """
    print(f"[*] Computing alpha diversity indices (Shannon, Simpson, Chao1) → {output_base}")

    # 1. Load and filter data
    df      = pd.read_excel(data_path, sheet_name=0)
    df_rank = df[df['Rank'].str.lower() == rank_level.lower()].copy()
    if df_rank.empty:
        raise ValueError(f"No data found for taxonomic level: {rank_level}")

    meta_cols   = {'Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name'}
    sample_cols = [c for c in df_rank.columns if c not in meta_cols]
    df_counts   = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

    # 2. Compute alpha diversity per sample from raw counts
    print(f"    [*] Computing diversity indices for {len(sample_cols)} sample(s)...")
    records = []
    for sample in sample_cols:
        metrics = compute_alpha_diversity(df_counts[sample])
        records.append({'Sample': sample, **metrics})

    alpha_df = pd.DataFrame(records)

    # 3. Map metadata groups
    meta_df = (pd.read_csv(metadata_path) if metadata_path.lower().endswith('.csv')
               else pd.read_excel(metadata_path))
    meta_df[sample_id_col] = meta_df[sample_id_col].astype(str).str.strip().str.lower()

    if pd.api.types.is_numeric_dtype(meta_df[category_col]):
        print(f"    [*] Continuous variable '{category_col}' — Quantile Binning (q=4).")
        meta_df[category_col] = pd.qcut(meta_df[category_col].dropna(), q=4).astype(str)

    meta_map = dict(zip(meta_df[sample_id_col], meta_df[category_col].astype(str)))

    alpha_df['CleanSample'] = alpha_df['Sample'].apply(
        lambda x: str(x).split('_')[0].strip().lower())
    alpha_df['Group'] = alpha_df['CleanSample'].map(meta_map)

    unmatched = alpha_df['Group'].isna().sum()
    if unmatched > 0:
        print(f"    [!] {unmatched} sample(s) not found in metadata — excluded.")
    alpha_df = alpha_df.dropna(subset=['Group'])

    if alpha_df.empty:
        raise ValueError("No sample IDs matched between abundance data and metadata.")

    # 4. Group order + color mapping from unified COLORS palette
    group_order = sorted(alpha_df['Group'].unique())
    palette_map = {g: COLORS[i % len(COLORS)] for i, g in enumerate(group_order)}

    # 5. Three-panel violin figure (Updated metric specs for stripped-down labels)
    metric_specs = [
        ('Shannon', "Shannon"),
        ('Simpson', 'Simpson'),
        ('Chao1',   'Chao1'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 7), facecolor='white')
    fig.patch.set_facecolor('white')

    for ax, (metric_key, metric_label) in zip(axes, metric_specs):
        data_by_group = [alpha_df[alpha_df['Group'] == g][metric_key].dropna().values
                         for g in group_order]

        # Violin bodies — only groups with > 1 observation
        valid_pos  = [i for i, d in enumerate(data_by_group) if len(d) > 1]
        valid_data = [data_by_group[i] for i in valid_pos]

        if valid_data:
            parts = ax.violinplot(valid_data, positions=valid_pos,
                                  widths=0.65, showmeans=False,
                                  showmedians=False, showextrema=False)
            for idx, body in zip(valid_pos, parts['bodies']):
                body.set_facecolor(palette_map[group_order[idx]])
                body.set_edgecolor('black')
                body.set_linewidth(1.1)
                body.set_alpha(0.60)

        # Mean ± SE point + error bar
        for i, g in enumerate(group_order):
            vals = data_by_group[i]
            if len(vals) == 0:
                continue
            mean_val = np.mean(vals)
            se_val   = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
            ax.errorbar(i, mean_val, yerr=se_val,
                        fmt='none', ecolor='black',
                        elinewidth=1.1, capsize=4, capthick=1.1, zorder=4)
            ax.plot(i, mean_val, 'o', color='black',
                    markersize=4.5, zorder=5,
                    markeredgecolor='black', markeredgewidth=0.7)

        # Removed title block as requested 
        # (Titles and KW p-values are no longer rendered on top of individual plots)

        # CLD letters (Kruskal-Wallis + Mann-Whitney Bonferroni)
        cld_input = alpha_df[['Group', metric_key]].rename(columns={metric_key: 'Value'})
        try:
            letters_dict = compute_kruskal_cld(cld_input, 'Value', 'Group')
        except Exception:
            letters_dict = {g: '' for g in group_order}

        y_vals = alpha_df[metric_key].dropna()
        y_max  = float(y_vals.max()) if len(y_vals) else 1.0
        y_pos  = y_max * 1.12
        ax.set_ylim(bottom=0, top=y_pos * 1.18)

        for i, g in enumerate(group_order):
            ax.text(i, y_pos, letters_dict.get(g, ''),
                    ha='center', va='bottom',
                    fontsize=11, fontweight='bold', color='black')

        # Panel styling — theme_bw
        ax.set_facecolor('white')
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)
            spine.set_color('black')
        ax.tick_params(axis='both', which='both', direction='out',
                       colors='black', width=0.8, length=4)
        ax.set_xticks(range(len(group_order)))
        ax.set_xticklabels(group_order, rotation=45, ha='right', fontsize=9)
        ax.set_xlabel(category_col, fontsize=11, fontweight='bold')
        ax.set_ylabel(metric_label,  fontsize=11, fontweight='bold')

    # 6. Shared legend (right side, moved nearer to the top via bbox_to_anchor adjustment)
    legend_patches = [
        mpatches.Patch(facecolor=palette_map[g], edgecolor='black',
                       linewidth=0.7, alpha=0.7, label=g)
        for g in group_order
    ]
    fig.legend(handles=legend_patches, title=category_col,
               title_fontsize=10, fontsize=10,
               loc='upper right', bbox_to_anchor=(1.0, 1.02),
               frameon=True, framealpha=0.95, edgecolor='#CCCCCC')

    plt.suptitle(f"{organism_name} — Alpha Diversity ({rank_level.capitalize()})",
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout(pad=2.0)
    export_figure(fig, output_base, fmt)
    plt.close(fig)
    print(f"    [✓] Alpha diversity plot saved: {output_base}.{fmt}")

    # 7. Summary table — Sample | Group | Shannon | Simpson | Chao1
    if not no_table:
        table_df = (
            alpha_df[['Sample', 'Group', 'Shannon', 'Simpson', 'Chao1']]
            .round({'Shannon': 4, 'Simpson': 4, 'Chao1': 2})
            .sort_values(['Group', 'Sample'])
            .reset_index(drop=True)
        )
        export_summary_table(table_df, output_base, sheet_name='Alpha_Diversity')
# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Violin Plot Generator — two modes:\n"
            "  taxa  : per-taxon relative abundance violins (ANOVA + Tukey HSD CLD)\n"
            "  alpha : Shannon H', Simpson (1-D), Chao1 violins (KW + Mann-Whitney CLD)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-d",    "--data",      required=True,
                        help="Path to taxa Excel file.")
    parser.add_argument("-r",    "--rank",      default="genus",
                        help="Taxonomic rank level.")
    parser.add_argument("-t",    "--threshold", type=float, default=0.01,
                        help="Mean relative abundance threshold for taxa mode (e.g. 0.01 = 1%%).")
    parser.add_argument("-o",    "--output",    default=None,
                        help="Output filename base (no extension).")
    parser.add_argument("-org",  "--organism",  default="Microbiome",
                        help="Prefix for the shared Y-axis label.")
    parser.add_argument("-m",    "--metadata",  required=True,
                        help="Path to metadata CSV or Excel file.")
    parser.add_argument("-c",    "--category",  required=True,
                        help="Metadata column to group by.")
    parser.add_argument("-id",   "--sample_id", default="SampleID",
                        help="Metadata column containing sample IDs.")
    parser.add_argument("-fmt",  "--format",    choices=["pdf", "png", "tiff"],
                        default="pdf", help="Output format.")
    parser.add_argument("--mode", choices=["taxa", "alpha"], default="taxa",
                        help=(
                            "Operational mode: "
                            "'taxa' = per-taxon relative abundance violin plots (ANOVA/Tukey); "
                            "'alpha' = Shannon, Simpson, Chao1 diversity index violins (KW/Mann-Whitney)."
                        ))
    parser.add_argument("--no_table", action="store_true",
                        help="Skip exporting the summary table (.xlsx).")

    args = parser.parse_args()

    OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'results', 'Violin_ANOVA')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    safe_name = args.organism.replace(" ", "_")

    if args.mode == "taxa":
        base_name   = (args.output if args.output
                       else f"{safe_name}_{args.rank}_{args.threshold}_Violin_ANOVA")
        output_base = os.path.join(OUTPUT_DIR, base_name)
        generate_taxa_violin_plots(
            data_path     = args.data,
            metadata_path = args.metadata,
            category_col  = args.category,
            sample_id_col = args.sample_id,
            rank_level    = args.rank,
            threshold     = args.threshold,
            organism_name = args.organism,
            output_base   = output_base,
            fmt           = args.format.lower(),
            no_table      = args.no_table,
        )

    elif args.mode == "alpha":
        base_name   = (args.output if args.output
                       else f"{safe_name}_{args.rank}_Alpha_Diversity")
        output_base = os.path.join(OUTPUT_DIR, base_name)
        generate_alpha_diversity_violins(
            data_path     = args.data,
            metadata_path = args.metadata,
            category_col  = args.category,
            sample_id_col = args.sample_id,
            rank_level    = args.rank,
            organism_name = args.organism,
            output_base   = output_base,
            fmt           = args.format.lower(),
            no_table      = args.no_table,
        )
