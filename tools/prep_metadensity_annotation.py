import argparse, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import gaussian_filter
from collections import defaultdict
import pyBigWig
import math
import pickle
import gzip
import os
import time

from metadensity_utils import (
    parse_gff3_for_boundaries,
    save_pickle
)

parser = argparse.ArgumentParser(description='Prep Metadensity Annotation - Extract gene boundaries from GFF3')

parser.add_argument('--GFF3', type=str, required=True,
                    help='Path to GFF3 annotation file (can be gzipped)')

args = parser.parse_args()

GFF3 = args.GFF3

print("=" * 80)
print("Metadensity Annotation Preparation")
print("=" * 80)
print(f"Input GFF3: {GFF3}")
print("=" * 80)

# Generate boundaries with full feature support
print("\nParsing GFF3 file...")
print("Extracting features:")
print("  - Exons (including first/last)")
print("  - Introns (computed from exons)")
print("  - 5' UTR and 3' UTR")
print("  - CDS (including first/last)")
print("\nThis may take several minutes for large annotation files...")

boundaries = parse_gff3_for_boundaries(GFF3)

print(f"\nParsing complete!")
print(f"  Processed {len(boundaries)} genes")

# Print feature statistics
feature_counts = defaultdict(int)
for gene_name, gene_features in boundaries.items():
    for feature_type in gene_features:
        feature_counts[feature_type] += 1

print("\nFeature statistics:")
for feature_type in sorted(feature_counts.keys()):
    print(f"  {feature_type}: {feature_counts[feature_type]} genes")

# Save boundaries
os.makedirs("./annotations", exist_ok=True)
output_path = os.path.join("./annotations/gene_boundaries.pkl")
save_pickle(boundaries, output_path)

print(f"\nSaved gene boundaries to: {output_path}")
print("=" * 80)