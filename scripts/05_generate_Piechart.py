import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.backends.backend_pdf import PdfPages

def generate_global_pie_chart(data_path, rank_level='domain', threshold=0.01, output_file='global_pie_chart.pdf'):
    """
    Generates a PDF pie chart of the global average relative abundance.
    Uses an external legend to prevent overlapping in highly skewed distributions.
    """
    # 1. Load data
    df = pd.read_excel(data_path, sheet_name=0)
    
    # 2. Filter strictly by the requested taxonomic rank
    df_rank = df[df['Rank'].str.lower() == rank_level.lower()].copy()
    if df_rank.empty:
        raise ValueError(f"No data found for the taxonomic level: {rank_level}")
        
    # Standardize the index
    if 'Scientific Name' in df_rank.columns:
        df_rank.set_index('Scientific Name', inplace=True)
    elif 'Name' in df_rank.columns:
        df_rank.set_index('Name', inplace=True)
    
    # 3. Isolate sample count columns
    exclude_cols = ['Rank', 'TaxID', 'original_header', 'Name', 'Scientific Name']
    sample_cols = [col for col in df_rank.columns if col not in exclude_cols]
    
    # Convert to pure numeric tensors
    df_counts = df_rank[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # 4. Mathematical Transformation: Relative Abundance per Sample
    rel_abund = df_counts.div(df_counts.sum(axis=0), axis=1)
    
    # 5. Calculate global arithmetic mean
    mean_rel_abund = rel_abund.mean(axis=1)
    
    # 6. Apply Thresholding
    high_abund = mean_rel_abund[mean_rel_abund >= threshold]
    low_abund = mean_rel_abund[mean_rel_abund < threshold]
    
    plot_data = high_abund.copy()
    
    # Aggregate low abundance taxa into "Others"
    if not low_abund.empty and low_abund.sum() > 0:
        plot_data['Others (<{:.1%})'.format(threshold)] = low_abund.sum()
        
    # Sort values descending for a monotonic visual gradient
    plot_data = plot_data.sort_values(ascending=False)
    
    # User's Custom Color Palette (30 Colors)
    colors = [
        "#b988d5", "#cbd588", "#88a05b", "#ffe156", "#6b8ba4", "#fda47d", 
        "#8e5634", "#d44645", "#5d9ca3", "#63b7af", "#dcd3ff", "#ff94cc", 
        "#ffa45b", "#806d40", "#2a363b", "#99b898", "#feceab", "#ff847c", 
        "#e84a5f", "#2a363b", "#56a5cc", "#c3e88d", "#ffcc5c", "#b09f59", 
        "#ff5e57", "#674d3c", "#4c4f69", "#8372a8", "#ff7c43", "#6a8a82"
    ]
    
    # Piecewise function to hide text inside slices smaller than 1.5% to avoid overlapping bounds
    def autopct_generator(pct):
        return f'{pct:.1f}%' if pct >= 1.5 else ''

    # 7. Render Pie Chart and export to PDF
    with PdfPages(output_file) as pdf:
        # Create a wider figure to comfortably fit the pie on the left and legend on the right
        fig, ax = plt.subplots(figsize=(14, 8))
        
        wedges, texts, autotexts = ax.pie(
            plot_data, 
            labels=None,  # Labels set to None. We move them to the legend.
            autopct=autopct_generator, 
            startangle=140,
            colors=colors[:len(plot_data)],
            wedgeprops=dict(edgecolor='white', linewidth=1.5)
        )
        
        # Style the percentage labels inside the pie
        plt.setp(autotexts, size=11, weight="bold", color="white")
        for autotext in autotexts:
            autotext.set_path_effects([path_effects.withStroke(linewidth=2, foreground='black')])
            
        # 8. Legend Engineering (Crucial step to preserve the < 1% data visually)
        # Create mathematically exact labels for the legend
        total_sum = plot_data.sum()
        legend_labels = [f"{idx} ({ (val/total_sum)*100:.2f}% )" for idx, val in plot_data.items()]
        
        # Position the legend completely outside the axes
        ax.legend(wedges, legend_labels,
                  title=f"{rank_level.capitalize()} Taxonomy",
                  title_fontproperties={'weight': 'bold', 'size': 14},
                  loc="center left",
                  bbox_to_anchor=(1.05, 0.5), # Anchor it to the right of the pie chart
                  fontsize=12,
                  frameon=False,
                  handlelength=1.5,
                  handleheight=1.5)
        
        ax.set_title(f"Global Relative Abundance ({rank_level.capitalize()})", 
                     fontsize=18, fontweight='bold', pad=20)
        
        # Ensure the bounding box captures the newly placed external legend
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight', pad_inches=0.5)
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Taxonomic Pie Chart Generator")
    parser.add_argument('-d', '--data', type=str, required=True, help='Path to taxa Excel file.')
    parser.add_argument('-r', '--rank', type=str, default='domain', help='Taxonomic rank (e.g., domain, phylum, species).')
    parser.add_argument('-t', '--threshold', type=float, default=0.01, help='Threshold to group taxa into "Others" (default: 0.01 for 1%).')
    parser.add_argument('-o', '--output', type=str, help='Custom output PDF filename.')

    args, unknown = parser.parse_known_args()
    
    out_file = args.output if args.output else f"PieChart_{args.rank}_{args.threshold}.pdf"
    print(f"Generating pie chart for rank '{args.rank}': {out_file}")
    
    generate_global_pie_chart(
        data_path=args.data,
        rank_level=args.rank,
        threshold=args.threshold,
        output_file=out_file
    )
    print("Pie chart successfully generated.")
