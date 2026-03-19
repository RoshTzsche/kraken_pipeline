import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.distance import pdist, squareform
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import scipy.stats as stats

COLORS = [

        "#b988d5", "#cbd588", "#88a05b", "#ffe156", "#6b8ba4", "#fda47d", 

        "#8e5634", "#d44645", "#5d9ca3", "#63b7af", "#dcd3ff", "#ff94cc", 

        "#ffa45b", "#806d40", "#2a363b", "#99b898", "#feceab", "#ff847c", 

        "#e84a5f", "#2a363b", "#56a5cc", "#c3e88d", "#ffcc5c", "#b09f59", 

        "#ff5e57", "#674d3c", "#4c4f69", "#8372a8", "#ff7c43", "#6a8a82"] 

def autopct_generator(pct):
    return f'{pct:.1f}%' if pct >= 1.5 else ''

def generate_global_pie_chart(df_rank, rank_level, threshold, output_file):
    """Phase 1: Generates the Relative Abundance Pie Chart."""
    print(f"[*] Generating Pie Chart -> {output_file}")
    
    sample_cols = [col for col in df_rank.columns if col not in ['Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name']]
    df_counts = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    rel_abund = df_counts.div(df_counts.sum(axis=0), axis=1)
    mean_rel_abund = rel_abund.mean(axis=1)
    
    high_abund = mean_rel_abund[mean_rel_abund >= threshold]
    low_abund = mean_rel_abund[mean_rel_abund < threshold]
    
    plot_data = high_abund.copy()
    if not low_abund.empty and low_abund.sum() > 0:
        plot_data['Others (<{:.1%})'.format(threshold)] = low_abund.sum()
        
    plot_data = plot_data.sort_values(ascending=False)
    
    with PdfPages(output_file) as pdf:
        fig, ax = plt.subplots(figsize=(14, 8))
        wedges, texts, autotexts = ax.pie(
            plot_data, labels=None, autopct=autopct_generator, 
            startangle=140, colors=COLORS[:len(plot_data)],
            wedgeprops=dict(edgecolor='white', linewidth=1.5)
        )
        
        plt.setp(autotexts, size=11, weight="bold", color="white")
        for autotext in autotexts:
            autotext.set_path_effects([path_effects.withStroke(linewidth=2, foreground='black')])
            
        total_sum = plot_data.sum()
        legend_labels = [f"{idx} ({ (val/total_sum)*100:.2f}% )" for idx, val in plot_data.items()]
        
        ax.legend(wedges, legend_labels, title=f"{rank_level.capitalize()} Taxonomy",
                  title_fontproperties={'weight': 'bold', 'size': 14}, loc="center left",
                  bbox_to_anchor=(1.05, 0.5), fontsize=12, frameon=False,
                  handlelength=1.5, handleheight=1.5)
        
        ax.set_title(f"Global Relative Abundance ({rank_level.capitalize()})", fontsize=18, fontweight='bold', pad=20)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight', pad_inches=0.5)
        plt.close()

def confidence_ellipse(x, y, ax, n_std=2.447, facecolor='none', **kwargs):
    """
    Mathematical Implementation of the Covariance Confidence Ellipse.
    n_std = 2.447 represents ~95% confidence for a 2D distribution (sqrt(chi2.ppf(0.95, 2)))
    """
    if x.size < 3:
        return  # Mathematically impossible to define a 2D area with < 3 points

    # Calculate covariance and pearson correlation
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    
    # Eigenvalues of the 2D correlation matrix
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Scale by standard deviation and translate to the empirical mean
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def generate_pcoa_plot(df_rank, rank_level, metadata_path, category_col, sample_id_col, output_file):
    """Phase 2: Generates the Principal Coordinate Analysis (PCoA) with Ellipses & Enhanced Design."""
    print(f"[*] Generating PCoA Plot -> {output_file}")
    
    # 1. Parse Data Matrix (Samples x Taxa)
    sample_cols = [col for col in df_rank.columns if col not in ['Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name']]
    df_counts = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    rel_abund = df_counts.div(df_counts.sum(axis=0), axis=1).T
    sample_names = rel_abund.index.tolist()
    
    # 2. Compute Bray-Curtis Distance Matrix
    dist_array = pdist(rel_abund.values, metric='braycurtis')
    dist_matrix = squareform(dist_array)
    
    # 3. Torgerson-Gower Scaling (Double Centering)
    n = dist_matrix.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n  
    B = -0.5 * J.dot(dist_matrix ** 2).dot(J)  
    
    # 4. Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    eigenvalues[eigenvalues < 0] = 0
    
    coords = eigenvectors[:, :2] * np.sqrt(eigenvalues[:2])
    variance_explained = (eigenvalues / np.sum(eigenvalues)) * 100
    pc1_var, pc2_var = variance_explained[0], variance_explained[1]
# 5. Robust Metadata Mapping & Continuous Data Binning
    groups = ["All Samples"] * n
    if metadata_path and category_col:
        meta_df = pd.read_csv(metadata_path) if metadata_path.endswith('.csv') else pd.read_excel(metadata_path)
        
        # --- AUTOMATIC NUMERICAL DISCRETIZATION (Binning) ---
        # If the column is numerical (like ORP, pH), bin it into 4 Quartiles
        if pd.api.types.is_numeric_dtype(meta_df[category_col]):
            print(f"    [*] Continuous numerical variable detected for '{category_col}'.")
            print(f"    [*] Applying Quantile Binning (q=4) to ensure sufficient group sizes for ellipses.")
            
            try:
                # pd.qcut partitions the data into 4 equal-sized probability bins
                bins = pd.qcut(meta_df[category_col].dropna(), q=4)
                meta_df[category_col] = bins.astype(str)
            except ValueError:
                # Fallback to equal-width binning if there are too many identical numbers
                bins = pd.cut(meta_df[category_col].dropna(), bins=4)
                meta_df[category_col] = bins.astype(str)
        # ----------------------------------------------------

        # Aggressively clean metadata IDs for robust matching
        meta_df[sample_id_col] = meta_df[sample_id_col].astype(str).str.strip().str.lower()
        
        # Create mapping dictionary, ensuring categories are cast as strings
        meta_map = dict(zip(meta_df[sample_id_col], meta_df[category_col].astype(str)))
        
        groups = []
        for name in sample_names:
            # Clean abundance sample name
            clean_name = name.split('_')[0].strip().lower()
            
            # Debug check to reveal the "Unknown" issue
            if clean_name not in meta_map:
                print(f"    [!] WARNING: Abundance Sample '{name}' (cleaned to '{clean_name}') NOT FOUND in metadata '{sample_id_col}' column.")
            
            # Map the value; if it evaluates to a string 'nan', label it Unknown
            mapped_val = meta_map.get(clean_name, "Unknown")
            if mapped_val == 'nan':
                mapped_val = "Unknown"
                
            groups.append(mapped_val)
            
    df_pcoa = pd.DataFrame({'PC1': coords[:, 0], 'PC2': coords[:, 1], 'Group': groups, 'Sample': sample_names})
    unique_groups = df_pcoa['Group'].unique()
    # 6. Render Enhanced PCoA Plot
    with PdfPages(output_file) as pdf:
        # Improved Design: Grey background standard for modern bioinformatics plots (like ggplot)
        fig, ax = plt.subplots(figsize=(11, 8), facecolor='#f8f9fa')
        ax.set_facecolor('#f4f4f6')
        ax.grid(color='white', linestyle='-', linewidth=1.5, alpha=0.8) # Beautiful clean grid
        
        for i, group in enumerate(unique_groups):
            subset = df_pcoa[df_pcoa['Group'] == group]
            color = COLORS[i % len(COLORS)]
            
            # Plot Points
            ax.scatter(subset['PC1'], subset['PC2'], 
                       s=180, alpha=0.9, label=group,
                       color=color, edgecolors='white', linewidth=2, zorder=3)
            
            # Plot 95% Confidence Ellipses (Requires at least 3 points per group)
            if len(subset) >= 3 and group != "Unknown":
                confidence_ellipse(subset['PC1'].values, subset['PC2'].values, ax, 
                                   n_std=2.447, edgecolor=color, facecolor=color, 
                                   alpha=0.15, linewidth=2, zorder=2)
                # Add dashed border to ellipse for style
                confidence_ellipse(subset['PC1'].values, subset['PC2'].values, ax, 
                                   n_std=2.447, edgecolor=color, facecolor='none', 
                                   linestyle='--', alpha=0.8, linewidth=1.5, zorder=2)

        # Draw axes origin lines
        ax.axhline(0, color='black', linestyle=':', linewidth=1.2, alpha=0.6, zorder=1)
        ax.axvline(0, color='black', linestyle=':', linewidth=1.2, alpha=0.6, zorder=1)
        
        # Typography & Styling
        ax.set_xlabel(f"PC1 ({pc1_var:.1f}%)", fontsize=14, fontweight='bold', color='#333333')
        ax.set_ylabel(f"PC2 ({pc2_var:.1f}%)", fontsize=14, fontweight='bold', color='#333333')
        ax.set_title(f"PCoA (Bray-Curtis) - {rank_level.capitalize()}", fontsize=18, fontweight='heavy', pad=20, color='#1a1a1a')
        
        # Remove borders for a cleaner look
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        # Legend styling
        leg = ax.legend(title=category_col if category_col else "Group", 
                        fontsize=12, title_fontproperties={'weight':'bold', 'size':13},
                        bbox_to_anchor=(1.03, 0.5), loc='center left', frameon=True,
                        facecolor='white', edgecolor='white', shadow=True)
        leg.get_title().set_color('#333333')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight', pad_inches=0.5, facecolor=fig.get_facecolor())
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Microbiome Visualization Suite: Pie Charts & PCoA")
    parser.add_argument('-d', '--data', type=str, required=True, help='Path to taxonomic_classification_clean.xlsx')
    parser.add_argument('-mode', '--mode', type=str, choices=['pie', 'pcoa', 'both'], default='both', help='Type of analysis to generate.')
    parser.add_argument('-r', '--rank', type=str, default='domain', help='Taxonomic rank (e.g., domain, genus).')
    
    # Pie Chart Specific
    parser.add_argument('-t', '--threshold', type=float, default=0.01, help='Pie Chart threshold (default: 0.01).')
    
    # PCoA Specific (Metadata)
    parser.add_argument('-m', '--metadata', type=str, help='Path to metadata CSV/Excel for PCoA coloring.')
    parser.add_argument('-c', '--category', type=str, help='Metadata column to group by in PCoA.')
    parser.add_argument('-id', '--sample_id', type=str, default='SampleID', help='Metadata column name for sample IDs.')

    args = parser.parse_args()
    
    # 1. Load the core dataset
    df = pd.read_excel(args.data, sheet_name=0)
    df_rank = df[df['Rank'].str.lower() == args.rank.lower()].copy()
    
    if df_rank.empty:
        raise ValueError(f"CRITICAL ERROR: No data found for the taxonomic rank: {args.rank}")
        
    if 'Scientific Name' in df_rank.columns:
        df_rank.set_index('Scientific Name', inplace=True)
    elif 'Name' in df_rank.columns:
        df_rank.set_index('Name', inplace=True)

    # 2. Execute Modes
    if args.mode in ['pie', 'both']:
        out_pie = f"PieChart_{args.rank}_{args.threshold}.pdf"
        generate_global_pie_chart(df_rank, args.rank, args.threshold, out_pie)
        
    if args.mode in ['pcoa', 'both']:
        # If metadata is used, append the category to the filename
        suffix = f"_{args.category}" if args.category else ""
        out_pcoa = f"PCoA_{args.rank}{suffix}.pdf"
        generate_pcoa_plot(df_rank, args.rank, args.metadata, args.category, args.sample_id, out_pcoa)
        
    print("[✓] Execution complete.")
