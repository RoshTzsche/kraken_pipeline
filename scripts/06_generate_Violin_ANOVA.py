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
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import math

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
# ggplot2-style discrete hue palette
# ---------------------------------------------------------------------------
def ggplot2_palette(n):
    import colorsys
    return [colorsys.hls_to_rgb(h / n, 0.65, 0.70) for h in range(n)]


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
def export_figure(fig, output_base, fmt, pad_inches=0.35):
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


# ---------------------------------------------------------------------------
# Compact Letter Display (Tukey HSD)
# ---------------------------------------------------------------------------
def calculate_cld(df, val_col, group_col):
    mean_series = df.groupby(group_col)[val_col].mean().sort_values(ascending=False)
    groups = mean_series.index.tolist()

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
# Draw one violin facet panel
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
# Main pipeline
# ---------------------------------------------------------------------------
def generate_taxa_violin_plots(data_path, metadata_path, category_col,
                                sample_id_col, rank_level, threshold,
                                organism_name, output_base, fmt):
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
    #    Right margin reserved for the vertical legend
    #    Top margin reserved for strip headers (per row)
    taxa_list  = list(melted["Taxon"].unique())
    n_vars     = len(taxa_list)
    ncols      = min(3, n_vars)
    nrows      = math.ceil(n_vars / ncols)

    legend_w_in = 1.8   # inches reserved on the right for the legend
    fig_w       = ncols * 4.8 + legend_w_in
    fig_h       = nrows * 5.2 + 0.7   # +0.7 bottom for shared x-label

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h),
                              facecolor="white")
    fig.patch.set_facecolor("white")

    if n_vars == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # 8. Draw each panel (no individual labels/titles)
    for i, taxon in enumerate(taxa_list):
        sub_df = melted[melted["Taxon"] == taxon].dropna(subset=["Abundance"])
        draw_violin_panel(axes[i], sub_df, group_order, palette_map)

    # Hide unused axes
    for j in range(n_vars, len(axes)):
        axes[j].set_visible(False)

    # 8b. X tick labels — only on the BOTTOM-MOST panel of each column.
    #     Handles incomplete last rows: every column independently finds its
    #     own lowest panel, so a 2-row column and a 3-row column each show
    #     labels only at their respective bottom.
    #
    #     col_panel_indices  →  list of flat indices that belong to this column,
    #                           ordered top → bottom (row 0, row 1, …)
    for col in range(ncols):
        col_panel_indices = [
            col + row * ncols
            for row in range(nrows)
            if col + row * ncols < n_vars        # skip out-of-range cells
        ]
        # All rows except the last: hide x tick labels but keep the spine/ticks
        for idx in col_panel_indices[:-1]:
            axes[idx].tick_params(axis="x", labelbottom=False)
            axes[idx].set_xticklabels([])        # belt-and-suspenders
        # Bottom row of this column: make sure labels are fully visible
        axes[col_panel_indices[-1]].tick_params(
            axis="x", labelbottom=True, labelrotation=45, labelsize=9
        )

    # 9. tight_layout with reserved margins
    #    rect = [left, bottom, right, top]
    legend_frac = legend_w_in / fig_w
    strip_frac  = (STRIP_H_IN * nrows + 0.3) / fig_h   # total vertical space for strips

    fig.tight_layout(
        rect=[0.04,                      # left  — supylabel room
              0.06,                      # bottom — supxlabel room
              1.0 - legend_frac - 0.01,  # right — legend room
              1.0 - strip_frac],         # top   — strip header room
        pad=1.6, h_pad=4.0, w_pad=2.0,
    )

    # 10. Shared axis labels — printed once per figure
    fig.supxlabel(category_col,
                  fontsize=12, fontweight="bold", y=0.01, va="bottom")
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
        # anchor at right edge of the plot area
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


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Violin Plot Generator (ANOVA + Tukey HSD) — R facet_wrap style"
    )
    parser.add_argument("-d",   "--data",      required=True,
                        help="Path to taxa Excel file.")
    parser.add_argument("-r",   "--rank",      default="genus",
                        help="Taxonomic rank level.")
    parser.add_argument("-t",   "--threshold", type=float, default=0.01,
                        help="Mean relative abundance threshold (e.g. 0.01 = 1%%).")
    parser.add_argument("-o",   "--output",    default=None,
                        help="Output filename base (no extension).")
    parser.add_argument("-org", "--organism",  default="Microbiome",
                        help="Prefix for the shared Y-axis label.")
    parser.add_argument("-m",   "--metadata",  required=True,
                        help="Path to metadata CSV or Excel file.")
    parser.add_argument("-c",   "--category",  required=True,
                        help="Metadata column to group by.")
    parser.add_argument("-id",  "--sample_id", default="SampleID",
                        help="Metadata column containing sample IDs.")
    parser.add_argument("-fmt", "--format",    choices=["pdf", "png", "tiff"],
                        default="pdf", help="Output format.")

    args = parser.parse_args()

    OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'Violin_ANOVA')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    safe_name   = args.organism.replace(" ", "_")
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
    )
