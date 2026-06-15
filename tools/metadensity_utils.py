from scipy.ndimage import gaussian_filter
import matplotlib.gridspec as gridspec
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pyBigWig
import pickle
import gzip
import math
import os

# Feature length for plots
feat_len = {
    'five_prime_UTR': 200,
    'three_prime_UTR': 300,
    'intron': 1000,
    'exon': 200,
    'first_exon': 200,
    'last_exon': 200,
    'CDS': 200,
    'first_CDS': 200,
    'last_CDS': 200,}

# Axis width ratios for plotting
ax_width_dict = {
    'UTR': 1,
    'exon': 2,
    'intron': 3,
    'CDS': 2,}

# Feature list
generic_rna = ['first_exon', 'exon', 'intron', 'last_exon']
protein_coding = ['five_prime_UTR', 'first_CDS', 'CDS', 'last_CDS', 'three_prime_UTR']

def _interval_overlap(a_start, a_end, b_start, b_end):
    # Returns overlap length in bp (0 if none). Assumes half-open [start, end).
    left = max(a_start, b_start)
    right = min(a_end, b_end)
    return max(0, right - left)

def load_flipper_results(filepath):
    df = pd.read_csv(filepath, sep='\t')
    return df

def _center_distance_to_interval(center, start, end):
    # 0 if inside; else distance to nearest edge.
    if start <= center <= end:
        return 0
    if center < start:
        return start - center
    return center - end

def filter_significant_windows(flipper_df, padj_threshold=0.05, l2fc_threshold=1.0):
    sig_df = flipper_df[flipper_df['padj'] <= padj_threshold].copy()
    upregulated = sig_df[sig_df['log2FoldChange'] >= l2fc_threshold].copy()
    downregulated = sig_df[sig_df['log2FoldChange'] <= -l2fc_threshold].copy()
    return upregulated, downregulated


def save_pickle(obj, filepath):
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'wb') as f:
        pickle.dump(obj, f)


def find_bigwig_files(bigwig_dir, control, treatment):
    """
    Find BigWig files for IP and INPUT samples

    Parameters:
    -----------
    bigwig_dir : str
        Directory containing BigWig files
    control : str
        Control sample name
    treatment : str
        Treatment sample name

    Returns:
    --------
    dict
        Dictionary with structure:
        {sample: {'IP': {'plus': path, 'minus': path},
                  'INPUT': {'plus': path, 'minus': path}}}
        or {sample: {'IP': path, 'INPUT': path}} for flat structure
    """
    bigwig_files = {}

    # Check if we have strand-split structure (plus/minus subdirectories)
    has_strand_dirs = (os.path.isdir(os.path.join(bigwig_dir, 'plus')) and
                       os.path.isdir(os.path.join(bigwig_dir, 'minus')))

    if has_strand_dirs:
        # Handle strand-split structure
        print(f"  Detected strand-split BigWig structure (plus/minus subdirectories)")

        for sample in [control, treatment]:
            bigwig_files[sample] = {}

            for condition in ['IP', 'INPUT']:
                bigwig_files[sample][condition] = {}

                for strand in ['plus', 'minus']:
                    strand_dir = os.path.join(bigwig_dir, strand)

                    # Set up patterns specific to each condition
                    # Pattern examples:
                    #   IP:    {sample}_IP_{replicate}.*.{strand}.bw
                    #   INPUT: {sample}_IN_{replicate}.*.{strand}.bw (abbreviated)
                    if condition == 'IP':
                        patterns = [
                            f"{sample}_IP_",
                            f"{sample}.IP."
                        ]
                    elif condition == 'INPUT':
                        patterns = [
                            f"{sample}_IN_", 
                            f"{sample}_INPUT_",
                            f"{sample}.INPUT.",
                            f"{sample}.IN."
                        ]

                    found = False
                    for fname in os.listdir(strand_dir):
                        if not fname.endswith('.bw'):
                            continue

                        # Check if filename matches any pattern and contains strand
                        fname_lower = fname.lower()
                        for pattern in patterns:
                            if pattern.lower() in fname_lower and strand in fname_lower:
                                bigwig_files[sample][condition][strand] = os.path.join(strand_dir, fname)
                                found = True
                                break

                        if found:
                            break

                # Validate that both strands were found
                if 'plus' not in bigwig_files[sample][condition] or 'minus' not in bigwig_files[sample][condition]:
                    raise FileNotFoundError(
                        f"Could not find both plus and minus strand BigWig files for {sample} {condition} in {bigwig_dir}\n"
                        f"Expected pattern: {sample}_{{IP/IN}}_{{rep}}.*.{{plus/minus}}.bw"
                    )

    else:
        # Handle flat structure (original behavior)
        print(f"  Detected flat BigWig structure")

        for sample in [control, treatment]:
            bigwig_files[sample] = {}

            # Find IP and INPUT files
            for condition in ['IP', 'INPUT']:
                # Look for files matching pattern: {sample}_{condition}.bw or similar
                pattern = f"{sample}_{condition}"
                for fname in os.listdir(bigwig_dir):
                    if pattern.lower() in fname.lower() and fname.endswith('.bw'):
                        bigwig_files[sample][condition] = os.path.join(bigwig_dir, fname)
                        break

            # Validate that both IP and INPUT were found
            if 'IP' not in bigwig_files[sample] or 'INPUT' not in bigwig_files[sample]:
                raise FileNotFoundError(
                    f"Could not find both IP and INPUT BigWig files for sample {sample} in {bigwig_dir}"
                )

    return bigwig_files


def parse_gff3_for_boundaries(gff3_file):
    """
    Parse GFF3 file to extract gene boundaries for all features

    Extracts:
    - Exons (with first/last designation)
    - Introns (derived from exons)
    - 5' UTR and 3' UTR
    - CDS (with first/last designation)

    Returns:
    --------
    dict
        Nested dictionary: {gene_name: {feature_type: [(start, end, strand), ...]}}
        Feature types: 'exon', 'first_exon', 'last_exon', 'intron',
                      'five_prime_UTR', 'three_prime_UTR', 'CDS', 'first_CDS', 'last_CDS'
    """
    boundaries = defaultdict(lambda: defaultdict(list))
    gene_exons = defaultdict(list)
    gene_info = {}

    open_func = gzip.open if gff3_file.endswith('.gz') else open

    with open_func(gff3_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)

            # Parse attributes to get gene_name
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            gene_name = attr_dict.get('gene_name') or attr_dict.get('gene_id')
            if not gene_name:
                continue

            # Store gene info
            if gene_name not in gene_info:
                gene_info[gene_name] = {'chrom': chrom, 'strand': strand}

            # Process different feature types
            if feature_type == 'gene':
                boundaries[gene_name]['gene'] = [(start, end, strand)]

            elif feature_type == 'exon':
                gene_exons[gene_name].append((start, end, strand))
                boundaries[gene_name]['exon'].append((start, end, strand))

            elif feature_type == 'five_prime_UTR' or feature_type == '5UTR':
                boundaries[gene_name]['five_prime_UTR'].append((start, end, strand))

            elif feature_type == 'three_prime_UTR' or feature_type == '3UTR':
                boundaries[gene_name]['three_prime_UTR'].append((start, end, strand))

            elif feature_type == 'CDS':
                boundaries[gene_name]['CDS'].append((start, end, strand))

    # derive introns, first/last exons, and first/last CDS
    for gene_name in gene_exons:
        exons = sorted(gene_exons[gene_name], key=lambda x: x[0])
        strand = exons[0][2]

        if len(exons) > 0:
            # Determine first and last exons based on strand
            if strand == '+':
                first_exon = exons[0]
                last_exon = exons[-1]
            else:  # minus strand
                first_exon = exons[-1]  # rightmost is first in transcription
                last_exon = exons[0]    # leftmost is last in transcription

            boundaries[gene_name]['first_exon'] = [first_exon]
            boundaries[gene_name]['last_exon'] = [last_exon]

        # Compute introns
        if len(exons) > 1:
            for i in range(len(exons) - 1):
                intron_start = exons[i][1] + 1  # End of current exon + 1
                intron_end = exons[i + 1][0] - 1  # Start of next exon - 1
                if intron_end > intron_start:  # Valid intron
                    boundaries[gene_name]['intron'].append((intron_start, intron_end, strand))

        # first/last CDS
        if 'CDS' in boundaries[gene_name] and len(boundaries[gene_name]['CDS']) > 0:
            cds_list = sorted(boundaries[gene_name]['CDS'], key=lambda x: x[0])

            if strand == '+':
                first_cds = cds_list[0]
                last_cds = cds_list[-1]
            else:  # minus strand
                first_cds = cds_list[-1]
                last_cds = cds_list[0]

            boundaries[gene_name]['first_CDS'] = [first_cds]
            boundaries[gene_name]['last_CDS'] = [last_cds]

    boundaries = {k: dict(v) for k, v in boundaries.items()}

    return boundaries

def annotate_flipper_windows_with_boundaries(windows_df, boundaries, max_distance=5000):
    annotated = defaultdict(list)

    for _, row in windows_df.iterrows():
        chrom = row['chr']          # Consider renaming / supporting 'chrom'
        win_start = int(row['start'])
        win_end = int(row['end'])
        gene_name = row['gene_name']

        if gene_name not in boundaries:
            continue

        gene_bounds = boundaries[gene_name]
        center = (win_start + win_end) // 2

        best = None
        # best = (overlap_bp, -dist_bp, feature_type, feat_start, feat_end, strand)
        # We’ll maximize overlap, and if tie minimize distance.

        for feature_type, feats in gene_bounds.items():
            if feature_type == 'gene':
                continue

            for feat_start, feat_end, strand in feats:
                # Quick reject by max_distance using interval expansion
                if (win_start > feat_end + max_distance) or (win_end < feat_start - max_distance):
                    continue

                overlap = _interval_overlap(win_start, win_end, feat_start, feat_end)
                dist = _center_distance_to_interval(center, feat_start, feat_end)

                if dist > max_distance:
                    continue

                candidate = (overlap, -dist, feature_type, feat_start, feat_end, strand)
                if (best is None) or (candidate > best):
                    best = candidate

        if best is None:
            continue

        overlap, neg_dist, feature_type, feat_start, feat_end, strand = best
        annotated[feature_type].append({
            'chrom': chrom,
            'start': win_start,
            'end': win_end,
            'gene_name': gene_name,
            'strand': strand,
            'feature_start': feat_start,
            'feature_end': feat_end,
            'window_center': center,
            'overlap_bp': overlap,
            'distance_bp': -neg_dist
        })

    return dict(annotated)

def dedup_annotated_by_feature(annotated_dict):
    """
    Deduplicate an annotated dict so each unique feature interval appears once.

    Keeps the 'best' entry per feature based on:
      - max overlap_bp
      - min distance_bp
      - max window length (end-start)
    """
    deduped = {}

    for feature_type, entries in annotated_dict.items():
        best_for_feature = {}

        for e in entries:
            key = (
                feature_type,
                e['chrom'],
                e['gene_name'],
                e['strand'],
                int(e['feature_start']),
                int(e['feature_end']),
            )

            overlap = int(e.get('overlap_bp', 0))
            dist = int(e.get('distance_bp', 10**12))
            win_len = int(e.get('end', 0)) - int(e.get('start', 0))

            score = (overlap, -dist, win_len)  # higher is better

            if key not in best_for_feature or score > best_for_feature[key][0]:
                best_for_feature[key] = (score, e)

        deduped[feature_type] = [v[1] for v in best_for_feature.values()]

    return deduped

def relative_information_transform(ip_signal, input_signal, truncate=False):
    raw = np.array(ip_signal, dtype=float)
    control_raw = np.array(input_signal, dtype=float)

    if truncate:
        ip_pseudocount = 0.0001
        control_pseudocount = 0.0001
    else:
        non_zero = raw[np.where(raw > 0)]
        if len(non_zero) == 0:
            ip_pseudocount = 1
        else:
            ip_pseudocount = np.min(non_zero)

        non_zero = control_raw[np.where(control_raw > 0)]
        if len(non_zero) == 0:
            control_pseudocount = 1 
        else:
            control_pseudocount = np.min(non_zero)

    raw_added = raw + ip_pseudocount
    control_added = control_raw + control_pseudocount
    ip_norm = raw_added / np.sum(raw_added)
    input_norm = control_added / np.sum(control_added)

    # Calculate relative information: pi * log2(pi/qi)
    values = np.multiply(
        ip_norm,
        np.log2(np.divide(ip_norm, input_norm))
    )
    return values


def trim_or_pad(density, target_length, align='left', pad_value=0):
    density = np.array(density)

    if len(density) == target_length:
        return density
    elif len(density) > target_length:
        # Trim
        if align == 'left':
            return density[:target_length]
        elif align == 'right':
            return density[-target_length:]
    else:
        # Pad
        discrepancy = target_length - len(density)
        if align == 'left':
            density = np.array(list(density) + [pad_value] * discrepancy)
        elif align == 'right':
            density = np.array([pad_value] * discrepancy + list(density))
        return density


def extract_flipper_densities_metadensity_style(
    bigwig_files,
    up_annotated,
    down_annotated,
    boundaries,
    window_size=100,
    features_to_extract=['exon', 'intron']):
    """
    Extract read densities from BigWig files and apply relative information transform

    Parameters:
    -----------
    bigwig_files : dict
        {sample_name: {'IP': path, 'INPUT': path}} for flat structure
        or {sample_name: {'IP': {'plus': path, 'minus': path}, 'INPUT': {...}}} for strand-split
    up_annotated : dict
        Upregulated windows by feature type
    down_annotated : dict
        Downregulated windows by feature type
    boundaries : dict
        Gene boundaries
    window_size : int
        Default window size for extraction
    features_to_extract : list
        List of features to extract densities for

    Returns:
    --------
    dict
        Nested dictionary: {regulation: {(feature, align, condition): np.array}}
        Shape of arrays: (n_windows, feature_length)
    """
    density_arrays = {
        'up': {},
        'down': {}
    }

    # Check if we have strand-split structure
    sample_name = list(bigwig_files.keys())[0]
    has_strand_split = isinstance(bigwig_files[sample_name]['IP'], dict)

    # Process each condition (up and down regulation)
    for regulation, annotated_dict in [('up', up_annotated), ('down', down_annotated)]:

        # Process each feature
        for feature in features_to_extract:
            if feature not in annotated_dict:
                continue

            windows = annotated_dict[feature]
            if len(windows) == 0:
                continue

            # Get feature length from config
            feature_length = feat_len.get(feature, window_size)

            # Initialize arrays for each alignment and condition
            for align in ['left', 'right']:
                for sample_name in bigwig_files.keys():

                    ip_densities = []
                    input_densities = []
                    ri_densities = []

                    # Open BigWig files based on structure
                    if has_strand_split:
                        # Strand-specific files - open both strands
                        ip_bw_plus = pyBigWig.open(bigwig_files[sample_name]['IP']['plus'])
                        ip_bw_minus = pyBigWig.open(bigwig_files[sample_name]['IP']['minus'])
                        input_bw_plus = pyBigWig.open(bigwig_files[sample_name]['INPUT']['plus'])
                        input_bw_minus = pyBigWig.open(bigwig_files[sample_name]['INPUT']['minus'])
                    else:
                        # Single file per condition
                        ip_bw = pyBigWig.open(bigwig_files[sample_name]['IP'])
                        input_bw = pyBigWig.open(bigwig_files[sample_name]['INPUT'])

                    # Extract density for each window
                    for window in windows:
                        chrom = window['chrom']
                        feat_start = window['feature_start']
                        feat_end = window['feature_end']
                        strand = window['strand']

                        try:
                            if has_strand_split:
                                # Use strand-specific file
                                if strand == '+':
                                    ip_signal = ip_bw_plus.values(chrom, feat_start, feat_end, numpy=True)
                                    input_signal = input_bw_plus.values(chrom, feat_start, feat_end, numpy=True)
                                else:  # minus strand
                                    ip_signal = ip_bw_minus.values(chrom, feat_start, feat_end, numpy=True)
                                    input_signal = input_bw_minus.values(chrom, feat_start, feat_end, numpy=True)
                                    # Reverse minus strand to 5'->3' orientation
                                    ip_signal = ip_signal[::-1]
                                    input_signal = input_signal[::-1]
                            else:
                                # Extract from combined file
                                ip_signal = ip_bw.values(chrom, feat_start, feat_end, numpy=True)
                                input_signal = input_bw.values(chrom, feat_start, feat_end, numpy=True)

                                # Reverse if on minus strand
                                if strand == '-':
                                    ip_signal = ip_signal[::-1]
                                    input_signal = input_signal[::-1]

                            # Handle NaN values
                            ip_signal = np.nan_to_num(ip_signal, nan=0.0)
                            input_signal = np.nan_to_num(input_signal, nan=0.0)

                            # trim or pad to align
                            ip_signal = trim_or_pad(ip_signal, feature_length, align=align)
                            input_signal = trim_or_pad(input_signal, feature_length, align=align)

                            # relative information
                            ri_signal = relative_information_transform(ip_signal, input_signal, truncate=False)

                            ip_densities.append(ip_signal)
                            input_densities.append(input_signal)
                            ri_densities.append(ri_signal)

                        except Exception as e:
                            # If extraction fails, add zeros
                            print(f"Warning: Failed to extract {chrom}:{feat_start}-{feat_end} ({strand}): {e}")
                            ip_densities.append(np.zeros(feature_length))
                            input_densities.append(np.zeros(feature_length))
                            ri_densities.append(np.zeros(feature_length))

                    if has_strand_split:
                        ip_bw_plus.close()
                        ip_bw_minus.close()
                        input_bw_plus.close()
                        input_bw_minus.close()
                    else:
                        ip_bw.close()
                        input_bw.close()

                    density_arrays[regulation][(feature, align, f'{sample_name}_IP')] = np.array(ip_densities)
                    density_arrays[regulation][(feature, align, f'{sample_name}_INPUT')] = np.array(input_densities)
                    density_arrays[regulation][(feature, align, f'{sample_name}_RI')] = np.array(ri_densities)

    return density_arrays



def calculate_grid_width(features_to_show, ax_width_dict):
    width = {}
    exist_keys = list(ax_width_dict.keys())

    for feat in features_to_show:
        valid_key = feat
        for key in exist_keys:
            if key in feat:
                valid_key = key
                break

        try:
            width[feat] = ax_width_dict[valid_key]
        except:
            width[feat] = ax_width_dict.get(feat, 2)

    return width


def get_xticklabels(feat, align, n_tick=5):
    flen = feat_len.get(feat, 100)

    if align == 'right':
        xticks = np.arange(0, flen + 1, flen / n_tick)
        xticklabel = ['{:.0f}'.format(x) for x in np.arange(-flen, 1, flen / 5)]

        # Special labels for 3' end
        if feat == 'three_prime_UTR':
            xticklabel[-1] = 'TTS'
        elif 'intron' in feat.split('_')[-1]:
            xticklabel[-1] = "3' SS"
        elif feat == 'last_CDS':
            xticklabel[-1] = 'stop codon'
    else:  # left alignment
        xticks = np.arange(0, flen, flen / 5)
        xticklabel = ['{:.0f}'.format(x) for x in np.arange(0, flen, flen / 5)]

        # Special labels for 5' end
        if feat == 'five_prime_UTR':
            xticklabel[0] = 'TSS'
        elif 'intron' in feat.split('_')[-1]:
            xticklabel[0] = "5' SS"
        elif feat == 'first_CDS':
            xticklabel[0] = 'start codon'

    return xticks, xticklabel


def generate_flipper_axis(nrows=2, features_to_show=None, height=None):
    if features_to_show is None:
        features_to_show = ['exon', 'intron']

    width_dict = calculate_grid_width(features_to_show, ax_width_dict)
    total_width = sum(list(width_dict.values())) * 2

    if height is None:
        height = nrows * 3

    fig = plt.figure(figsize=(total_width, height))
    spec = gridspec.GridSpec(ncols=total_width, nrows=nrows, figure=fig)

    ax_dict = {}

    for r in range(nrows):
        current = 0
        for feat in features_to_show:
            for align in ['left', 'right']:
                width = width_dict[feat]

                # Create subplot with shared y-axis
                if feat == features_to_show[0] and align == 'left':
                    ax_dict[feat, align, 'rep{}'.format(r + 1)] = fig.add_subplot(
                        spec[r, current:current + width]
                    )
                else:
                    ax_dict[feat, align, 'rep{}'.format(r + 1)] = fig.add_subplot(
                        spec[r, current:current + width],
                        sharey=ax_dict[features_to_show[0], 'left', 'rep{}'.format(r + 1)]
                    )
                    plt.setp(ax_dict[feat, align, 'rep{}'.format(r + 1)].get_yticklabels(),
                            visible=False)

                # Add title to first row
                if r == 0:
                    lbl = "5'" if align == 'left' else "3'"
                    featstr = feat.replace('five_prime_UTR', "5' UTR").replace('three_prime_UTR', "3' UTR")
                    ax_dict[feat, align, 'rep{}'.format(r + 1)].set_title(f'{featstr}({lbl})')

                # Hide xticklabels for non-last rows
                if r + 1 < nrows:
                    plt.setp(ax_dict[feat, align, 'rep{}'.format(r + 1)].get_xticklabels(),
                            visible=False)

                # Set x-tick labels
                xticks, xticklabel = get_xticklabels(feat, align, n_tick=5)
                ax_dict[feat, align, 'rep{}'.format(r + 1)].set_xticks(xticks)
                ax_dict[feat, align, 'rep{}'.format(r + 1)].set_xticklabels(xticklabel, rotation=90)

                current = current + width

    return fig, ax_dict



def plot_metadensity(
    density_arrays,
    up_counts,
    down_counts,
    features_to_show=['exon', 'intron'],
    output_path=None,
    smooth=False,
    sigma=5,
    alpha=0.3,
    ymax=None,
    plot_type='mean_density',
    conditions_to_plot=None,
    font_scale=1.4,
    legend_labels=None,
    legend_titles="Samples"):
    if plot_type == 'mean_density':
        return plot_mean_density_flipper(
            density_arrays, up_counts, down_counts,
            features_to_show, smooth, sigma, alpha, ymax, conditions_to_plot, font_scale, legend_labels, legend_titles)
    elif plot_type == 'heatmap':
        return plot_heatmap_flipper(
            density_arrays, features_to_show, conditions_to_plot)
    else:
        raise ValueError(f"Unknown plot_type: {plot_type}")

def plot_mean_density_flipper(
    density_arrays,
    up_counts,
    down_counts,
    features_to_show=['exon', 'intron'],
    smooth=False,
    sigma=5,
    alpha=0.3,
    ymax=None,
    conditions_to_plot=None,
    font_scale=1.4,
    legend_labels=None,
    legend_title="Samples"):

    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.ndimage import gaussian_filter

    # -------------------------
    # Global font scaling
    # -------------------------
    base_font = 10 * font_scale
    plt.rcParams.update({
        "font.size": base_font,
        "axes.titlesize": base_font * 1.2,
        "axes.labelsize": base_font * 1.1,
        "xtick.labelsize": base_font,
        "ytick.labelsize": base_font,
        "legend.fontsize": base_font,
        "legend.title_fontsize": base_font * 1.1
    })

    # Generate axes
    fig, ax_dict = generate_flipper_axis(nrows=2, features_to_show=features_to_show)

    # Define color scheme
    color_dict = {
        'up': '#d62728',
        'down': '#1f77b4'
    }

    if conditions_to_plot is None:
        sample_keys = set()
        for regulation in density_arrays.keys():
            for key in density_arrays[regulation].keys():
                feature, align, condition = key
                sample_keys.add(condition)
        conditions_to_plot = list(sample_keys)

    row_info = {}

    # -------------------------
    # Plot up/down
    # -------------------------
    for row_idx, regulation in enumerate(['up', 'down']):

        if regulation not in density_arrays:
            continue

        den_arr = density_arrays[regulation]

        for feat in features_to_show:
            for align in ['left', 'right']:

                ax = ax_dict[feat, align, f'rep{row_idx + 1}']

                for condition in conditions_to_plot:

                    key = (feat, align, condition)
                    if key not in den_arr:
                        continue

                    density_data = den_arr[key]
                    if density_data.shape[0] == 0:
                        continue

                    md = np.nanmean(density_data, axis=0)
                    std = np.nanstd(density_data, axis=0)
                    n = density_data.shape[0]
                    sem = std / np.sqrt(n)

                    if smooth:
                        md = gaussian_filter(md, sigma=sigma)
                        sem = gaussian_filter(sem, sigma=sigma)

                    # Custom legend label mapping
                    label = legend_labels.get(condition, condition) if legend_labels else condition

                    ax.plot(md, label=label, linewidth=2)
                    ax.fill_between(
                        np.arange(len(md)),
                        md - sem,
                        md + sem,
                        alpha=alpha)

                if ymax is not None:
                    ax.set_ylim(ymin=0, ymax=ymax)

                if feat == features_to_show[0] and align == 'left':
                    ylabel = 'Relative Information' if any('RI' in c for c in conditions_to_plot) else 'Read Density'
                    ax.set_ylabel(ylabel)

        total_windows = sum(up_counts.get(feat, 0) for feat in features_to_show) if regulation == 'up' else sum(down_counts.get(feat, 0) for feat in features_to_show)
        row_info[regulation] = total_windows

    # -------------------------
    # Legend
    # -------------------------
    handles, labels = ax_dict[features_to_show[0], 'left', 'rep1'].get_legend_handles_labels()

    legend = fig.legend(
        handles,
        labels,
        bbox_to_anchor=(0.85, 0.7),
        loc='upper left',
        ncol=1,
        title=legend_title,
        frameon=True,
        fancybox=True)

    # -------------------------
    # Info box
    # -------------------------
    info_text = []
    if 'up' in row_info:
        info_text.append(f"Row 1: Upregulated\n(n={row_info['up']} windows)")
    if 'down' in row_info:
        info_text.append(f"Row 2: Downregulated\n(n={row_info['down']} windows)")

    if info_text:
        info_box_text = '\n\n'.join(info_text)
        fig.text(
            0.86,
            0.35,
            info_box_text,
            transform=fig.transFigure,
            fontsize=base_font,
            verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8, edgecolor='gray'))

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    return fig

def plot_heatmap_flipper(
    density_arrays,
    features_to_show=['exon', 'intron'],
    conditions_to_plot=None,
    cmap='RdBu_r',
    vmin=None,
    vmax=None):

    fig, ax_dict = generate_flipper_axis(nrows=2, features_to_show=features_to_show)
    if conditions_to_plot is None:
        sample_keys = set()
        for regulation in density_arrays.keys():
            for key in density_arrays[regulation].keys():
                feature, align, condition = key
                sample_keys.add(condition)
        conditions_to_plot = [list(sample_keys)[0]]

    condition = conditions_to_plot[0]

    # Plot for up / down
    for row_idx, regulation in enumerate(['up', 'down']):

        if regulation not in density_arrays:
            continue

        den_arr = density_arrays[regulation]

        # Plot each feature
        for feat in features_to_show:
            for align in ['left', 'right']:

                ax = ax_dict[feat, align, f'rep{row_idx + 1}']

                key = (feat, align, condition)
                if key not in den_arr:
                    continue

                density_data = den_arr[key]

                if density_data.shape[0] == 0:
                    continue

                im = ax.imshow(
                    density_data,
                    aspect='auto',
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    interpolation='none')

                if feat == features_to_show[0] and align == 'left':
                    ax.set_ylabel('Windows')

    plt.colorbar(im, ax=list(ax_dict.values()), label='Signal')
    plt.tight_layout()
    return fig