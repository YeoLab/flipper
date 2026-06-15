import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import gaussian_filter
from collections import defaultdict###
import pickle
import pyBigWig

from metadensity_utils import (
    load_flipper_results,
    filter_significant_windows,
    parse_gff3_for_boundaries,
    annotate_flipper_windows_with_boundaries,
    save_pickle,
    find_bigwig_files,
    generate_flipper_axis,
    extract_flipper_densities_metadensity_style,
    plot_metadensity,
    dedup_annotated_by_feature
)

parser = argparse.ArgumentParser(description='Window Metadensity with Metadensity-style visualization')

parser.add_argument('--diff_window_result', type=str, required=True,
                    help='Path to flipper differential window results (TSV)')
parser.add_argument('--diff_gene_result', type=str, required=True,
                    help='Path to flipper differential gene results (TSV)')
parser.add_argument('--boundaries_file', type=str, required=True,
                    help='Path to gene boundaries pickle file')
parser.add_argument('--control', type=str, required=True,
                    help='Control sample name')
parser.add_argument('--treatment', type=str, required=True,
                    help='Treatment sample name')
parser.add_argument('--BIGWIG_DIR', type=str, required=True,
                    help='Directory containing BigWig files')

parser.add_argument('--log2fc_threshold', type=float, default=1.0,
                    help='Log2 fold change threshold (default: 1.0)')
parser.add_argument('--padj_threshold', type=float, default=0.05,
                    help='Adjusted p-value threshold (default: 0.05)')
parser.add_argument('--top_n_genes', type=int, default=100,
                    help='Number of top genes to analyze by fisher p-value (default: 100)')

parser.add_argument('--features', type=str, default=None,
                    help='Comma-separated list of features (e.g., "exon,intron,five_prime_UTR")')
parser.add_argument('--smooth', action='store_true',
                    help='Apply Gaussian smoothing to density plots')
parser.add_argument('--sigma', type=float, default=5,
                    help='Smoothing sigma parameter (default: 5)')
parser.add_argument('--alpha', type=float, default=0.2,
                    help='Transparency for SEM bands (default: 0.2)')
parser.add_argument('--ymax', type=float, default=None,
                    help='Maximum y-axis value for plots')

parser.add_argument('--conditions', type=str, default='RI',
                    help='Comma-separated conditions to plot: IP, INPUT, RI (default: RI)')
parser.add_argument('--window_size', type=int, default=100,
                    help='Window size for density extraction (default: 100)')

args = parser.parse_args()

diff_window_result = args.diff_window_result
diff_gene_result = args.diff_gene_result
boundaries_file = args.boundaries_file
control = args.control
treatment = args.treatment
BIGWIG_DIR = args.BIGWIG_DIR
log2fc_threshold = args.log2fc_threshold
padj_threshold = args.padj_threshold
top_n_genes = args.top_n_genes
features_to_show = args.features.split(",")

conditions_to_plot = []
for cond in args.conditions.split(','):
    cond = cond.strip()
    if cond == 'RI':
        conditions_to_plot.extend([f'{control}_RI', f'{treatment}_RI'])
    elif cond == 'IP':
        conditions_to_plot.extend([f'{control}_IP', f'{treatment}_IP'])
    elif cond == 'INPUT':
        conditions_to_plot.extend([f'{control}_INPUT', f'{treatment}_INPUT'])

print("=" * 80)
print("flipper Metadensity Analysis Pipeline")
print("=" * 80)
print(f"Input files:")
print(f"  - flipper windows: {diff_window_result}")
print(f"  - flipper genes: {diff_gene_result}")
print(f"  - Boundaries: {boundaries_file}")
print(f"  - BigWig directory: {BIGWIG_DIR}")
print(f"\nSamples:")
print(f"  - Control: {control}")
print(f"  - Treatment: {treatment}")
print(f"\nFilters:")
print(f"  - Top N genes: {top_n_genes}")
print(f"  - padj threshold: {padj_threshold}")
print(f"  - log2FC threshold: {log2fc_threshold}")
print(f"\nFeatures to analyze: {', '.join(features_to_show)}")
print(f"Conditions to plot: {', '.join(conditions_to_plot)}")
print("=" * 80)

print("\n[1/5] Loading flipper results...")
flipper_df = load_flipper_results(diff_window_result)
print(f"  Loaded {len(flipper_df)} differential windows")

upregulated = flipper_df[flipper_df['log2FoldChange'] >= 0].copy()
downregulated = flipper_df[flipper_df['log2FoldChange'] <= 0].copy()

print(f"  Upregulated: {len(upregulated)} windows")
print(f"  Downregulated: {len(downregulated)} windows")

print("\n[2/5] Loading gene boundaries and annotating windows...")
with open(boundaries_file, 'rb') as f:
    boundaries = pickle.load(f)
print(f"  Loaded boundaries for {len(boundaries)} genes")

# Locate the features that contain/are closest too each of the windows tested. 
up_annotated = annotate_flipper_windows_with_boundaries(upregulated, boundaries, max_distance=5000)
down_annotated = annotate_flipper_windows_with_boundaries(downregulated, boundaries, max_distance=5000)

# Remove duplciates (sometimes windows are right next to each other, these features should not be counted twice). 
up_annotated = dedup_annotated_by_feature(up_annotated)
down_annotated = dedup_annotated_by_feature(down_annotated)

save_pickle(up_annotated, "./annotations/up_regulate_window_annotated.pkl")
save_pickle(down_annotated, "./annotations/down_regulate_window_annotated.pkl")
print("\n  Saved annotations to:", "./annotations/")

print("\n[3/5] Finding BigWig files...")
bigwig_files = find_bigwig_files(BIGWIG_DIR, control, treatment)
for sample, files in bigwig_files.items():
    print(f"  {sample}:")
    print(f"    IP: {files['IP']}")
    print(f"    INPUT: {files['INPUT']}")

print("\n[4/5] Extracting densities from BigWig files...")
print("  This may take several minutes depending on the number of windows...")
density_arrays = extract_flipper_densities_metadensity_style(
    bigwig_files,
    up_annotated,
    down_annotated,
    boundaries,
    window_size=args.window_size,
    features_to_extract=features_to_show
)
print("  Density extraction complete!")

save_pickle(density_arrays, "./annotations/window_densities.pkl")
print("  Saved density arrays to:", "./annotations/window_densities.pkl")

print("\n[5/5] Generating plots...")
up_counts = {k: len(v) for k, v in up_annotated.items()}
down_counts = {k: len(v) for k, v in down_annotated.items()}

os.makedirs("./plots", exist_ok=True)

fig = plot_metadensity(
    density_arrays,
    up_counts,
    down_counts,
    features_to_show=features_to_show,
    smooth=5,
    sigma=0.2,
    alpha=0.2,
    ymax=None,
    conditions_to_plot=conditions_to_plot,
    font_scale=1.4,
    legend_labels=None,
    legend_titles="Samples"
)
output_path = "./plots/metadensity_mean.svg"
fig.savefig(output_path, format="svg", bbox_inches="tight")
plt.show()
print(f"  Saved mean density plot to: {output_path}")