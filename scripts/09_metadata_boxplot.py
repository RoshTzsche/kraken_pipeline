import argparse
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

def get_cld_letters(df, val_col, group_col):
    """Calculates Tukey's HSD Compact Letter Display (CLD) for grouping significance."""
    # Group by treatment and calculate mean to sort groups from highest to lowest
    mean_series = df.groupby(group_col)[val_col].mean().sort_values(ascending=False)
    groups = mean_series.index.tolist()

    # Extract data arrays per group
    group_data = [df[df[group_col] == g][val_col].dropna().values for g in groups]
    
    # Check if there's enough data to perform statistics
    if len(groups) < 2 or df[val_col].nunique() <= 1:
        return {g: "a" for g in groups}

    # One-way ANOVA to check for global variance differences
    _, p_val = stats.f_oneway(*group_data)
    if p_val > 0.05 or pd.isna(p_val):
        return {g: "a" for g in groups} # No significant global difference

    # Tukey HSD post-hoc test
    tukey = pairwise_tukeyhsd(endog=df[val_col], groups=df[group_col], alpha=0.05)
    
    # Create a significance matrix based on Tukey results
    sig_matrix = pd.DataFrame(False, index=groups, columns=groups)
    for row in tukey.summary().data[1:]:
        g1, g2, reject = row[0], row[1], row[-1]
        if str(reject).strip().lower() == "true":
            sig_matrix.loc[g1, g2] = True
            sig_matrix.loc[g2, g1] = True

    # Assign non-overlapping letters to non-significantly different groups
    letters = {g: "" for g in groups}
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

    return letters

def main():
    parser = argparse.ArgumentParser(description="Generate a normalized grouped Boxplot with ANOVA and Tukey HSD from Metadata")
    parser.add_argument("-m", "--metadata", required=True, help="Path to metadata file (Excel or CSV)")
    parser.add_argument("-x", "--treatment", required=True, help="X-axis group column name (e.g. Type)")
    parser.add_argument("-y", "--morphology", required=True, nargs='+', help="List of numerical columns to plot together")
    parser.add_argument("--normalize", choices=["yes", "no"], default="yes",
                        help="Z-score normalize columns before plotting (default: yes). "
                             "Use 'no' to plot raw values with original units.")
    parser.add_argument("--no_legend", action="store_true",
                        help="Hide the color legend from the plot.")
    args = parser.parse_args()

    print(f"[*] Loading metadata from: {args.metadata}")
    
    # Read metadata file based on extension
    if args.metadata.endswith('.csv'):
        df = pd.read_csv(args.metadata)
    else:
        df = pd.read_excel(args.metadata)

    # Ensure all target columns exist in the dataframe
    missing_cols = [col for col in [args.treatment] + args.morphology if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in metadata: {missing_cols}")

    # Drop rows missing treatment data and enforce string type for grouping
    df = df.dropna(subset=[args.treatment])
    df[args.treatment] = df[args.treatment].astype(str)

    # Clean and normalize numerical columns
    normalize = args.normalize == "yes"
    if normalize:
        print("[*] Normalizing morphology columns (Z-score scaling)...")
    else:
        print("[*] Skipping normalization — plotting raw values...")
    plot_cols = []
    for col in args.morphology:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        if normalize:
            norm_col_name = f"{col}_Zscore"
            df[norm_col_name] = (df[col] - df[col].mean()) / df[col].std()
            plot_cols.append(norm_col_name)
        else:
            plot_cols.append(col)
    normalized_cols = plot_cols

    # Melt DataFrame to long format for Seaborn grouped boxplot compatibility
    df_melted = df.melt(id_vars=[args.treatment], 
                        value_vars=normalized_cols, 
                        var_name='Morphology_Trait', 
                        value_name='Normalized_Value')
    
    # Drop NaNs generated during melting
    df_melted = df_melted.dropna(subset=['Normalized_Value'])

    # Prepare data arrays for plotting iteration
    treatments = sorted(df_melted[args.treatment].unique())
    traits = normalized_cols
    n_treatments = len(treatments)

    # Calculate Tukey letters per morphology trait independently
    letters_dict = {}
    for trait in traits:
        trait_df = df_melted[df_melted['Morphology_Trait'] == trait]
        letters_dict[trait] = get_cld_letters(trait_df, 'Normalized_Value', args.treatment)

    # Plot configuration (publication style)
    sns.set_theme(style="ticks")
    dodge_width = 0.8

    if normalize:
        # --- NORMALIZED MODE: single figure, all traits on a shared Z-score Y-axis ---
        fig, ax_single = plt.subplots(figsize=(12, 6))

        try:
            ax = sns.boxplot(
                data=df_melted, x='Morphology_Trait', y='Normalized_Value',
                hue=args.treatment, hue_order=treatments,
                palette="Set2", width=dodge_width, gap=0.15,
                linewidth=1.5, flierprops=dict(marker='o', markersize=5, alpha=0.5),
                ax=ax_single
            )
        except TypeError:
            ax = sns.boxplot(
                data=df_melted, x='Morphology_Trait', y='Normalized_Value',
                hue=args.treatment, hue_order=treatments,
                palette="Set2", width=0.6, linewidth=1.5,
                flierprops=dict(marker='o', markersize=5, alpha=0.5),
                ax=ax_single
            )

        offsets = np.linspace(-dodge_width/2 + dodge_width/(2*n_treatments),
                               dodge_width/2 - dodge_width/(2*n_treatments),
                               n_treatments)
        y_max_global  = df_melted['Normalized_Value'].max()
        y_range        = y_max_global - df_melted['Normalized_Value'].min()
        text_y_offset  = y_range * 0.05

        for i, trait in enumerate(traits):
            trait_df = df_melted[df_melted['Morphology_Trait'] == trait]
            for j, treatment in enumerate(treatments):
                group_data = trait_df[trait_df[args.treatment] == treatment]['Normalized_Value']
                if not group_data.empty:
                    x_pos = i + offsets[j]
                    y_pos = group_data.max() + text_y_offset
                    ax.text(x_pos, y_pos, letters_dict[trait].get(treatment, ""),
                            ha='center', va='bottom', fontsize=10, fontweight='bold', color='black')

        ax.set_ylim(bottom=df_melted['Normalized_Value'].min() - text_y_offset,
                    top=y_max_global + (text_y_offset * 3))
        ax.set_xticklabels([col.replace('_Zscore', '') for col in traits], fontweight='bold')
        plt.xlabel("", fontsize=12, fontweight='bold')
        plt.ylabel("Z-score", fontsize=12, fontweight='bold')
        if not args.no_legend:
            plt.legend(title=args.treatment, title_fontsize='11', fontsize='10',
                       bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
        else:
            plt.legend().remove()
        sns.despine()
        plt.tight_layout()

    else:
        # --- RAW MODE: one subplot per variable, each with its own Y-axis scale ---
        n_traits  = len(traits)
        fig, axes = plt.subplots(1, n_traits, figsize=(5 * n_traits, 6), sharey=False)
        if n_traits == 1:
            axes = [axes]

        for idx, (trait, col_name) in enumerate(zip(traits, args.morphology)):
            ax = axes[idx]
            trait_df = df_melted[df_melted['Morphology_Trait'] == trait].copy()

            try:
                sns.boxplot(
                    data=trait_df, x=args.treatment, y='Normalized_Value',
                    hue=args.treatment, legend=False,
                    order=treatments, palette="Set2", width=dodge_width, gap=0.15,
                    linewidth=1.5, flierprops=dict(marker='o', markersize=5, alpha=0.5),
                    ax=ax
                )
            except TypeError:
                sns.boxplot(
                    data=trait_df, x=args.treatment, y='Normalized_Value',
                    hue=args.treatment, legend=False,
                    order=treatments, palette="Set2", width=0.6, linewidth=1.5,
                    flierprops=dict(marker='o', markersize=5, alpha=0.5),
                    ax=ax
                )

            # CLD letters — each subplot has its own Y scale
            y_max_local   = trait_df['Normalized_Value'].max()
            y_range_local = y_max_local - trait_df['Normalized_Value'].min()
            text_y_offset = y_range_local * 0.05

            for j, treatment in enumerate(treatments):
                group_data = trait_df[trait_df[args.treatment] == treatment]['Normalized_Value']
                if not group_data.empty:
                    ax.text(j, group_data.max() + text_y_offset,
                            letters_dict[trait].get(treatment, ""),
                            ha='center', va='bottom', fontsize=10, fontweight='bold', color='black')

            ax.set_ylim(bottom=trait_df['Normalized_Value'].min() - text_y_offset,
                        top=y_max_local + (text_y_offset * 3))
            ax.set_title(col_name, fontsize=12, fontweight='bold')
            ax.set_xlabel("")
            ax.set_ylabel(col_name, fontsize=11, fontweight='bold')
            ax.set_xticks(range(len(treatments)))
            ax.set_xticklabels(treatments, rotation=45, ha='right', fontweight='bold')

            sns.despine(ax=ax)

        # Single shared X-axis label across all subplots
        fig.supxlabel(args.treatment, fontsize=12, fontweight='bold')

        # Original legend style on the last subplot
        if not args.no_legend:
            plt.legend(title=args.treatment, title_fontsize='11', fontsize='10',
                       bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
        else:
            plt.legend().remove()

        plt.tight_layout()

    # 3. Dynamic output folder and filename generation
    # Extract only the first word/acronym of each morphology column to keep filenames readable
    safe_names = [col.split(' ')[0].replace('/', '_').replace('\\', '_') for col in args.morphology]
    combined_name = "_".join(safe_names)
    
    # Create output directory automatically inside results/
    output_dir = os.path.join("..", "results", f"Morphology_{combined_name}")
    os.makedirs(output_dir, exist_ok=True)
    
    output_filename = f"{combined_name}_normalized_boxplot.png"
    output_path = os.path.join(output_dir, output_filename)

    # Save to disk
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[✓] Plot successfully saved to: {output_path}")

if __name__ == "__main__":
    main()