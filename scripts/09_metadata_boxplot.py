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
    print("[*] Normalizing morphology columns (Z-score scaling)...")
    normalized_cols = []
    for col in args.morphology:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Z-score normalization: (x - mean) / std
        norm_col_name = f"{col}_Zscore"
        df[norm_col_name] = (df[col] - df[col].mean()) / df[col].std()
        normalized_cols.append(norm_col_name)

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
    plt.figure(figsize=(12, 6))

    # 1. Draw Grouped Boxplot
    dodge_width = 0.8
    try:
        # For Seaborn >= 0.13.0 (supports 'gap' parameter to separate hue boxes)
        ax = sns.boxplot(
            data=df_melted, 
            x='Morphology_Trait', 
            y='Normalized_Value', 
            hue=args.treatment,
            hue_order=treatments,
            palette="Set2", 
            width=dodge_width,
            gap=0.15, # <--- Adds separation between boxes of the same category
            linewidth=1.5,
            flierprops=dict(marker='o', markersize=5, alpha=0.5)
        )
    except TypeError:
        # Fallback for older Seaborn versions (reduces width to create space)
        ax = sns.boxplot(
            data=df_melted, 
            x='Morphology_Trait', 
            y='Normalized_Value', 
            hue=args.treatment,
            hue_order=treatments,
            palette="Set2", 
            width=0.6, 
            linewidth=1.5,
            flierprops=dict(marker='o', markersize=5, alpha=0.5)
        )

    # 2. Add significance letters above each corresponding box
    # Calculate offset positions for dodged boxes (Seaborn default logic)
    offsets = np.linspace(-dodge_width/2 + dodge_width/(2*n_treatments), 
                           dodge_width/2 - dodge_width/(2*n_treatments), 
                           n_treatments)

    # Determine a uniform Y offset for text placement based on global data range
    y_max_global = df_melted['Normalized_Value'].max()
    y_range = y_max_global - df_melted['Normalized_Value'].min()
    text_y_offset = y_range * 0.05

    # Iterate through traits and treatments to place CLD text
    for i, trait in enumerate(traits):
        trait_df = df_melted[df_melted['Morphology_Trait'] == trait]
        for j, treatment in enumerate(treatments):
            group_data = trait_df[trait_df[args.treatment] == treatment]['Normalized_Value']
            if not group_data.empty:
                y_max_group = group_data.max()
                letra = letters_dict[trait].get(treatment, "")
                
                # Position calculation: base x-tick (i) + dodge offset (offsets[j])
                x_pos = i + offsets[j]
                y_pos = y_max_group + text_y_offset
                
                ax.text(x_pos, y_pos, letra, ha='center', va='bottom', 
                        fontsize=10, fontweight='bold', color='black')

    # Formatting axes limits to prevent cut-off text
    ax.set_ylim(bottom=df_melted['Normalized_Value'].min() - text_y_offset, 
                top=y_max_global + (text_y_offset * 3))

    # Clean trait names for the X-axis labels (remove the "_Zscore" suffix)
    clean_labels = [col.replace('_Zscore', '') for col in traits]
    ax.set_xticklabels(clean_labels, fontweight='bold')
    
    # Setup Minimalist Labels (Removed Title and X-axis label as requested)
    plt.xlabel("", fontsize=12, fontweight='bold')
    plt.ylabel("Z-score", fontsize=12, fontweight='bold')
    
    # Legend formatting outside the plot
    plt.legend(title=args.treatment, title_fontsize='11', fontsize='10', 
               bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)

    sns.despine() # Clean borders
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