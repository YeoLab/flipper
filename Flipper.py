############################ Setup ##################################

# Import necessary packages. 
import subprocess, logging, time, re, sys, argparse, os, glob, json
from itertools import combinations
import pandas as pd, numpy as np
from functools import reduce
from time import sleep
import yaml

# Load up the main Flipper config file. 
locals().update(config)

# Load up the config file from the Skipper run.  
skipper_config = config["SKIPPER_CONFIG"]

if skipper_config.endswith(".py"):
    include: skipper_config  
elif skipper_config.endswith(".json") or skipper_config.endswith(".job"):
    with open(skipper_config) as f:
        cfg = json.load(f)
    for k, v in cfg.items():
        globals()[k] = v
elif skipper_config.endswith(".yaml"):
    with open(skipper_config) as f:
        cfg = yaml.safe_load(f)
        locals().update(cfg)
else:
    raise ValueError(f"Unsupported Skipper config file type: {skipper_config}")

############################ Process config file ##################################

# Load in inputs from the config. 
skipper_output = config["SKIPPER_OUTPUT"]
if skipper_output[-1] != "/":
    skipper_output = skipper_output + "/"
WORKDIR = config.get("WORKDIR")
workdir: config['WORKDIR']
control = config["CONTROL"]
treatment = config["TREATMENT"]

# Load configs with defaults. 
normalization = config.get("NORMALIZATION", "MOR_hier")
substrate_control = config.get("SUBSTRATE_CONTROL", "INf")
gc_method = config.get("GC_METHOD", "NONE")
length_method = config.get("LENGTH_METHOD", "NONE")
between_method = config.get("BETWEEN_METHOD", "NONE")
make_bigwigs = config.get("MAKE_BIGWIGS", True)
homer = config.get("HOMER", False)
metadensity = config.get("METADENSITY", False)
features = config.get("FEATURES", "first_exon,exon,intron,last_exon")
pval_threshold = config.get("PVAL_THRESHOLD", 0.05)
log2fc_threshold = config.get("LOG2FC_THRESHOLD", 1)

# Complete some minor processing from the Skipper config files. 
CHROM_SIZES = STAR_DIR + "/chrNameLength.txt"
UNINFORMATIVE_READ = str(3 - INFORMATIVE_READ)

# Create log file and tmp file paths. 
if not os.path.exists("stderr"): os.makedirs("stderr")
if not os.path.exists("stdout"): os.makedirs("stdout")
if not os.path.exists("tmp"): os.makedirs("tmp")

############################ Create window type to bam dictionaries ##################################
# (This is required for the homer/finemapping section). 

# Define sample labels.
conditions = [control, treatment]

# Define window types and link to samples (for homer)
window_type = ["downreg", "upreg"]
condition_to_windows = dict(zip(conditions, window_type))

# Load manifest, clean whitespace
manifest = (
    pd.read_csv(MANIFEST, comment="#")
    .dropna(subset=["Experiment", "Sample"])
    .apply(lambda col: col.str.strip() if col.dtype == "object" else col)
)

# Subset to just the control and treatment groups. 
manifest = manifest[manifest["Sample"].isin(conditions)]

# Build replicate labels
manifest["Input_replicate_label"] = manifest["Sample"].astype(str) + "_IN_" + manifest["Input_replicate"].astype(int).astype(str)
manifest["CLIP_replicate_label"]  = manifest["Sample"].astype(str) + "_IP_" + manifest["CLIP_replicate"].astype(int).astype(str)

# Use provided bams if present, otherwise get path to skipper generated bams. 
if {"Input_bam", "CLIP_bam"}.issubset(manifest.columns):
    manifest["Input_bam"] = manifest["Input_bam"].str.strip()
    manifest["CLIP_bam"]  = manifest["CLIP_bam"].str.strip()
else:
    manifest["Input_bam"] = skipper_output + "secondary_results/bams/dedup/genome/" + manifest["Input_replicate_label"] + ".genome.Aligned.sort.dedup.bam"
    manifest["CLIP_bam"]  = skipper_output + "secondary_results/bams/dedup/genome/" + manifest["CLIP_replicate_label"] + ".genome.Aligned.sort.dedup.bam"

# Collapse to dicts keyed by sample
sample_to_experiment = dict(zip(manifest["Sample"], manifest["Experiment"]))
sample_to_input_bams = manifest.groupby("Sample")["Input_bam"].apply(list).to_dict()
sample_to_clip_bams  = manifest.groupby("Sample")["CLIP_bam"].apply(list).to_dict()

# Map window type -> input BAMs.
window_to_input_bams = {
    window: sample_to_input_bams[exp]
    for exp, window in condition_to_windows.items()
}

# Map window type -> clip BAMs.
window_to_clip_bams = {
    window: sample_to_clip_bams[exp]
    for exp, window in condition_to_windows.items()
}

############################ Create helper functions ##################################

# Adjusts naming conventions to be consistent with R's X-prefix behavior. 
def _canon_sample_name(s: str) -> str:
    s = str(s)
    s = re.sub(r"^X(?=\d)", "", s)  
    s = s.replace(".", "-")        
    return s

# Used to extract scaling.
def get_div_factor_from_table(wc, input):
    path = input.scale_table

    df = pd.read_csv(path, sep="\t", comment="#")

    # Required columns.
    if "Sample" not in df.columns or "SizeFactor" not in df.columns:
        raise KeyError(
            f"Scale table at {path} must contain columns 'Sample' and 'SizeFactor'. "
            f"Found columns: {list(df.columns)}"
        )

    # Canonicalize both table samples and wildcard sample to handle R's X-prefix behavior.
    df["Sample_canon"] = df["Sample"].astype(str).map(_canon_sample_name)
    query_sample = _canon_sample_name(wc.sample)

    row = df.loc[df["Sample_canon"] == query_sample, "SizeFactor"]
    if row.empty:
        examples = ", ".join(df["Sample_canon"].dropna().unique()[:10])
        raise ValueError(
            f"No SizeFactor found for sample={wc.sample} (canon={query_sample}) in {path}. "
            f"Example samples: {examples}"
        )

    f = float(row.iloc[0])
    if f <= 0:
        raise ValueError(f"SizeFactor must be > 0 for sample={wc.sample}. Got {f}")

    # scaled = unscaled / SizeFactor
    return str(f)

############################ Creating rule all ##################################

homer_outputs = []
if homer:
    homer_outputs = [
        expand("intermediates/for_homer/finemapping/finemapping{window_type}.nt_coverage.bed",
               window_type=window_type),
        expand("intermediates/for_homer/region_matched_background/fixed/{window_type}.sampled_fixed_windows.bed.gz",
               window_type=window_type),
        expand("intermediates/for_homer/finemapping/mapped_sites/{window_type}.finemapped_windows.bed.gz",
               window_type=window_type),
        expand("homer/finemapped_results/{window_type}/homerResults.html",
               window_type=window_type),
    ]
    

metadensity_outputs = []
if metadensity:
    metadensity_outputs = [
        "plots/metadensity_mean.svg",
        "annotations/up_regulate_window_annotated.pkl",
        "annotations/down_regulate_window_annotated.pkl",
    ]

sample_labels = (
    pd.concat([
        manifest["Input_replicate_label"],
        manifest["CLIP_replicate_label"],
    ], ignore_index=True)
    .astype(str)
    .tolist()
)

# Normalized bedgraphs are only possible with MOR normalization. 
bedgraph_outputs = []
if normalization == "MOR_hier" and make_bigwigs:
    bedgraph_outputs = (
        expand("tracks/bedgraphs/{sample}.plus.scaled.bg",  sample=sample_labels) +
        expand("tracks/bedgraphs/{sample}.minus.scaled.bg", sample=sample_labels) +
        expand("tracks/bigwigs/{sample}.plus.scaled.bw",  sample=sample_labels) +
        expand("tracks/bigwigs/{sample}.minus.scaled.bw", sample=sample_labels)
    )

# Run rule all. 
rule all:
    input:
        window_counts_INsub = "intermediates/window_counts_INsub.tsv",
        norm_factors = "intermediates/normalization_factors.tsv",
        all_results = "tables/all_windows.tsv",
        diff_INsub_result = "tables/differential_windows.tsv",
        diff_gene_result = "tables/differential_genes.tsv",
        PCA = "plots/PCA.svg",
        feature = "plots/FeatureBar.svg",
        volcano = "plots/Volcano.png",
        MA = "plots/MAplot.png",
        Disp = "plots/DispEst.png",
        GO_plot = "plots/Gene_ontology.svg",
        *homer_outputs,
        *metadensity_outputs,
        *bedgraph_outputs

################################ Running sub-rules ####################################

rule preproccess:
    input:
        count_tables=expand(skipper_output + "secondary_results/counts/genome/tables/{experiment_label}.tsv.gz",experiment_label=manifest["Experiment"])
    output:
        window_counts_INsub = "intermediates/window_counts_INsub.tsv",
        all_count_data = "intermediates/all_count_data.tsv"
    params:
        control_label = sample_to_experiment[control],
        treatment_label = sample_to_experiment[treatment],
        skipper_output = skipper_output,
        anno_path = ancient(FEATURE_ANNOTATIONS),
        substrate_control = substrate_control,
        tool_dir = config['TOOLS']
    threads: 2
    resources:
        mem_mb = 16000,
        runtime = 30,
    log:
        stdout = "stdout/preproccess.out",
        stderr = "stderr/preproccess.err",
    benchmark: "benchmarks/preproccess/preproccess.txt"
    conda:
        "envs/skipper_R.yaml"
    shell:
        """
        Rscript --vanilla {params.tool_dir}/preprocess.R \
            {params.control_label} \
            {params.treatment_label} \
            {params.skipper_output} \
            {params.anno_path} \
            {params.substrate_control} \
            {params.tool_dir} \
        > {log.stdout} 2> {log.stderr}
        """

rule calc_norm_factors:
    input:
        count_tables=expand(skipper_output + "secondary_results/counts/genome/tables/{experiment_label}.tsv.gz",experiment_label=manifest["Experiment"]),
        all_count_data = "intermediates/all_count_data.tsv"
    output:
        norm_factors = "intermediates/normalization_factors.tsv",
        gene_wise_IN = "intermediates/gene_wise_IN.tsv",
    params:
        method = normalization, 
        control_label = sample_to_experiment[control],
        treatment_label = sample_to_experiment[treatment],
        skipper_output = skipper_output,
        anno_path = ancient(FEATURE_ANNOTATIONS),
        tool_dir = config['TOOLS'],
        gc_method = gc_method,
        length_method = length_method,
        between_method = between_method,
        substrate_control = substrate_control,
    threads: 1 
    resources:
        mem_mb = 16000,
        runtime = 30,
    log:
        stdout = "stdout/calc_norm_factors.out",
        stderr = "stderr/calc_norm_factors.err",
    benchmark: "benchmarks/preproccess/calc_norm.txt"
    conda:
        "envs/Flipper_main.yaml"
    shell:
        """
        Rscript --vanilla {params.tool_dir}/calc_norm_factors.R \
            {params.method} \
            {params.control_label} \
            {params.treatment_label} \
            {params.skipper_output} \
            {params.anno_path} \
            {params.tool_dir} \
            {params.gc_method} \
            {params.length_method} \
            {params.between_method} \
            {params.substrate_control} \
        > {log.stdout} 2> {log.stderr}
        """
         
rule differential_analysis:
    input:
        window_counts_INsub = "intermediates/window_counts_INsub.tsv",
        norm_factors = "intermediates/normalization_factors.tsv"
    output:
        all_result = "tables/all_windows.tsv",
        diff_INsub_result = "tables/differential_windows.tsv",
        diff_gene_result = "tables/differential_genes.tsv",
        INsub_PCA = "plots/PCA.svg",
        INsub_volcano = "plots/Volcano.png",
        INsub_feature = "plots/FeatureBar.svg",
        MA = "plots/MAplot.png",
        Disp = "plots/DispEst.png",
        windows = expand("intermediates/{window_type}_windows.tsv",window_type=window_type),
        all_result_norm = "tables/all_windows_norm_counts.tsv",
        diff_INsub_result_norm = "tables/differential_windows_norm_counts.tsv"
    params:
        norm_method = normalization, 
        pval_threshold = pval_threshold,
        log2fc_threshold = log2fc_threshold,
        substrate_control = substrate_control,
        tool_dir = config['TOOLS']
    threads: 2 
    resources:
        mem_mb = 16000,
        runtime = 30,
    log:
        stdout = "stdout/differential_analysis.out",
        stderr = "stderr/differential_analysis.err",
    benchmark: "benchmarks/preproccess/differential_analysis.txt"
    conda:
        "envs/Flipper_main.yaml"
    shell:
        """
        Rscript --vanilla {params.tool_dir}/differential_analysis.R \
            {input.window_counts_INsub} \
            {input.norm_factors} \
            {params.norm_method} \
            {params.pval_threshold} \
            {params.log2fc_threshold} \
            {params.substrate_control} \
            {params.tool_dir} \
        > {log.stdout} 2> {log.stderr}
        """

# This could potentially reuse Skipper nt_coverage by subsetting counts,
# but finemapping depends on regenerated window_n groups, so we currently
# recount from BAMs, just to be safe.
rule get_nt_coverage:
    input:
        windows    = "intermediates/{window_type}_windows.tsv",
        clip_bams  = lambda wildcards: window_to_clip_bams[wildcards.window_type],
        input_bams = lambda wildcards: window_to_input_bams[wildcards.window_type]
    output:
        nt_census       = temp("intermediates/for_homer/finemapping/finemapping{window_type}.nt_census.bed"),
        nt_input_counts = temp("intermediates/for_homer/finemapping/finemapping{window_type}.nt_coverage.input.counts"),
        nt_clip_counts  = temp("intermediates/for_homer/finemapping/finemapping{window_type}.nt_coverage.clip.counts"),
        nt_coverage     = "intermediates/for_homer/finemapping/finemapping{window_type}.nt_coverage.bed",
    threads: 2 
    resources:
        mem_mb = 48000,
        runtime = 60,
    log:
        stdout = "stdout/{window_type}.get_nt_coverage.out",
        stderr = "stderr/{window_type}.get_nt_coverage.err",
    benchmark: "benchmarks/get_nt_coverage/{window_type}.all_replicates.reproducible.txt"
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail
    
        # Count rows (skip header).
        nrows=$( (zcat -f "{input.windows}" | tail -n +2 | wc -l) || echo 0 )
    
        if [ "$nrows" -eq 0 ]; then
            echo "[get_nt_coverage] {wildcards.window_type}: no sites, writing empty outputs."
            mkdir -p "$(dirname {output.nt_census})"
            : > {output.nt_census}
            : > {output.nt_input_counts}
            : > {output.nt_clip_counts}
            : > {output.nt_coverage}
        else
            echo "Running on node: $(hostname)"
            echo "[`date`] Starting get_nt_coverage for {wildcards.window_type}"
    
            # Build per-nt census.
            zcat -f "{input.windows}" \
                | tail -n +2 \
                | sort -k1,1 -k2,2n -T tmp \
                | awk -v OFS="\t" '{{start = $2-37; if (start < 0) {{start = 0}}; print $1, start, $3+37, $4, $5, $6}}' \
                | bedtools merge -i - -s -c 6 -o distinct \
                | awk -v OFS="\t" '{{for (i = $2; i < $3; i++) {{print $1, i, i+1, "MW:" NR ":" i-$2, 0, $4, NR}} }}' \
                > {output.nt_census}
    
            # Input coverage counts.
            samtools cat {input.input_bams} \
                | bedtools intersect -s -wa -a - -b {output.nt_census} \
                | bedtools bamtobed -i - \
                | awk '($1 != "chrEBV") && ($4 !~ "/{UNINFORMATIVE_READ}$")' \
                | bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - \
                | bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - \
                | bedtools sort -i - \
                | bedtools coverage -counts -s -a {output.nt_census} -b - \
                | awk '{{print $NF}}' \
                > {output.nt_input_counts}
    
            # CLIP coverage counts.
            samtools cat {input.clip_bams} \
                | bedtools intersect -s -wa -a - -b {output.nt_census} \
                | bedtools bamtobed -i - \
                | awk '($1 != "chrEBV") && ($4 !~ "/{UNINFORMATIVE_READ}$")' \
                | bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - \
                | bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - \
                | bedtools sort -i - \
                | bedtools coverage -counts -s -a {output.nt_census} -b - \
                | awk '{{print $NF}}' \
                > {output.nt_clip_counts}
    
            # Final combined coverage table.
            paste {output.nt_census} {output.nt_input_counts} {output.nt_clip_counts} \
                > {output.nt_coverage}
    
            echo "[`date`] Finished get_nt_coverage for {wildcards.window_type}"
        fi \
        > {log.stdout} 2> {log.stderr}
        """

rule sample_background_windows_by_region:
    input:
        windows     = "intermediates/{window_type}_windows.tsv",
        all_windows = ancient(FEATURE_ANNOTATIONS),
    output:
        variable_windows = "intermediates/for_homer/region_matched_background/variable/{window_type}.sampled_variable_windows.bed.gz",
        fixed_windows    = "intermediates/for_homer/region_matched_background/fixed/{window_type}.sampled_fixed_windows.bed.gz"
    params:
        output_dir = "intermediates/for_homer/region_matched_background",
    threads: 1 
    resources:
        mem_mb = 16000,
        runtime = 10,
    log:
        stdout = "stdout/{window_type}.sample_background_windows_by_region.out",
        stderr = "stderr/{window_type}.sample_background_windows_by_region.err",
    benchmark: "benchmarks/sample_background_windows_by_region/{window_type}.sample_background_windows_by_region.txt"
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail
        nrows=$( (tail -n +2 "{input.windows}" | wc -l) || echo 0 )
        if [ "$nrows" -eq 0 ]; then
            echo "[sample_background_windows_by_region] {wildcards.window_type}: no sites, writing empty outputs."
            mkdir -p $(dirname {output.variable_windows}) $(dirname {output.fixed_windows})
            : | gzip -c > {output.variable_windows}
            : | gzip -c > {output.fixed_windows}
        else
            echo "[sample_background_windows_by_region] {wildcards.window_type}: found $nrows sites, running R script."
            Rscript --vanilla {TOOL_DIR}/sample_matched_background_by_region.R \
                {input.windows} {input.all_windows} 75 {params.output_dir} {wildcards.window_type}
        fi \
        > {log.stdout} 2> {log.stderr}
        """

rule finemap_windows:
    input:
        nt_coverage = "intermediates/for_homer/finemapping/finemapping{window_type}.nt_coverage.bed",
    output:
        finemapped_windows = "intermediates/for_homer/finemapping/mapped_sites/{window_type}.finemapped_windows.bed.gz"
    params:
        output_dir = "intermediates/for_homer/finemapping/mapped_sites",
    threads: 2 
    resources:
        mem_mb = 32000,
        runtime = 60,
    log:
        stdout = "stdout/{window_type}.finemap_windows.out",
        stderr = "stderr/{window_type}.finemap_windows.err",
    benchmark: "benchmarks/finemap_windows/{window_type}.all_replicates.reproducible.txt"
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail
        nrows=$( (grep -v '^#' "{input.nt_coverage}" | wc -l) || echo 0 )
        if [ "$nrows" -eq 0 ]; then
            echo "[finemap_windows] {wildcards.window_type}: no coverage rows, writing empty output."
            mkdir -p $(dirname {output.finemapped_windows})
            : | gzip -c > {output.finemapped_windows}
        else
            echo "[finemap_windows] {wildcards.window_type}: found $nrows rows, running R script."
            Rscript --vanilla {TOOL_DIR}/finemap_enriched_windows.R \
                {input.nt_coverage} {params.output_dir} {wildcards.window_type}
        fi \
        > {log.stdout} 2> {log.stderr}
        """

rule run_homer:
    input:
        finemapped_windows = "intermediates/for_homer/finemapping/mapped_sites/{window_type}.finemapped_windows.bed.gz",
        background         = "intermediates/for_homer/region_matched_background/fixed/{window_type}.sampled_fixed_windows.bed.gz",
        genome             = GENOME
    output:
        report = "homer/finemapped_results/{window_type}/homerResults.html"
    params:
        output_dir   = "homer/finemapped_results/{window_type}",
        preparse_dir = "homer/preparsed",
    threads: 4
    resources:
        mem_mb = 2000,
        runtime = 40,
    log:
        stdout = "stdout/{window_type}.run_homer.out",
        stderr = "stderr/{window_type}.run_homer.err",
    benchmark: "benchmarks/run_homer/{window_type}.all_replicates.reproducible.txt"
    conda:
        "envs/homer.yaml"
    shell:
        r"""
        set -euo pipefail
        nrows=$( (zcat -f "{input.finemapped_windows}" 2>/dev/null | grep -v '^#' | wc -l | tr -d '[:space:]') || echo 0 )

        if [ "$nrows" -eq 0 ]; then
            echo "[run_homer] {wildcards.window_type}: no finemapped windows, writing placeholder HTML."
            mkdir -p $(dirname {output.report})
            cat > {output.report} <<'HTML'
<!DOCTYPE html>
<meta charset="utf-8">
<title>No sites</title>
<p>No sites available; HOMER was skipped.</p>
HTML
        else
            echo "[run_homer] {wildcards.window_type}: found $nrows windows, running HOMER."
            findMotifsGenome.pl \
              <(zcat {input.finemapped_windows} | awk -F'\t' -v OFS='\t' '{{print $4 ":" $9, $1, $2+1, $3, $6}}') \
              {input.genome} {params.output_dir} \
              -preparsedDir {params.preparse_dir} -size given -rna -nofacts -S 20 -len 5,6,7,8,9 -nlen 1 \
              -bg <(zcat {input.background} | awk -F'\t' -v OFS='\t' '{{print $4, $1, $2+1, $3, $6}}')
        fi \
        > {log.stdout} 2> {log.stderr}
        """

rule Gene_ontology:
    input:
        diff_result = "tables/differential_genes.tsv",
        gene_wise_IN = "intermediates/gene_wise_IN.tsv"
    output:
        GO_plot = "plots/Gene_ontology.svg",
    params:
        species = config['SPECIES'],
        tool_dir = config['TOOLS']
    threads: 1
    resources:
        mem_mb = 12000,
        runtime = 30,
    log:
        stdout = "stdout/Gene_ontology.out",
        stderr = "stderr/Gene_ontology.err",
    benchmark: "benchmarks/postprocess/Gene_ontology.txt",
    conda: 
        "envs/R_GO.yaml"
    shell:
        """
        Rscript --vanilla {params.tool_dir}/gene_ontology.R \
            "{params.species}" \
            "{input.diff_result}" \
            "{input.gene_wise_IN}"\
            "{output.GO_plot}" \
            "{params.tool_dir}" \
        > {log.stdout} 2> {log.stderr}
        """
    
rule Prep_Metadensity_Annotation:
    input:
        GFF3 = GFF,
    output:
        gene_boundaries = "annotations/gene_boundaries.pkl",
    params:
        log2fc_threshold = log2fc_threshold,
        padj_threshold = pval_threshold,
        tool_dir = config['TOOLS']
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 30,
    log:
        stdout = "stdout/Prep_Metadensity_Annotation.out",
        stderr = "stderr/Prep_Metadensity_Annotation.err",
    benchmark: "benchmarks/postprocess/Metadensity.txt",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {params.tool_dir}/prep_metadensity_annotation.py \
            --GFF3 {input.GFF3} \
        > {log.stdout} 2> {log.stderr}
        """
        
rule Window_Metadensity:
    input:
        diff_window_result = "tables/differential_windows.tsv",
        diff_gene_result = "tables/differential_genes.tsv",
        boundaries_file = "annotations/gene_boundaries.pkl",
    output:
        Metadensity = "plots/metadensity_mean.svg",
        up_regulate_window_annotated = "annotations/up_regulate_window_annotated.pkl",
        down_regulate_window_annotated = "annotations/down_regulate_window_annotated.pkl",
    params:
        BIGWIG_DIR = config['SKIPPER_OUTPUT'] + "/secondary_results/bigwigs/scaled",
        control = config['CONTROL'],
        treatment = config['TREATMENT'],
        features = features, 
        log2fc_threshold = log2fc_threshold,
        padj_threshold = pval_threshold,
        tool_dir = config['TOOLS']
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 60,
    log:
        stdout = "stdout/Metadensity.out",
        stderr = "stderr/Metadensity.err",
    benchmark: "benchmarks/postprocess/Metadensity.txt",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {params.tool_dir}/metadensity.py \
            --diff_window_result {input.diff_window_result} \
            --diff_gene_result {input.diff_gene_result} \
            --boundaries_file {input.boundaries_file} \
            --control {params.control} \
            --treatment {params.treatment} \
            --BIGWIG_DIR {params.BIGWIG_DIR} \
            --log2fc_threshold {params.log2fc_threshold} \
            --padj_threshold {params.padj_threshold} \
            --features {params.features} \
        > {log.stdout} 2> {log.stderr}
        """

rule make_scaled_bigwig_from_bedgraph:
    input:
        bg_plus = skipper_output + "secondary_results/bedgraphs/unscaled/plus/{sample}.unscaled.plus.bg",
        bg_minus = skipper_output + "secondary_results/bedgraphs/unscaled/minus/{sample}.unscaled.minus.bg",
        chrom_sizes = STAR_DIR + "/chrNameLength.txt",
        scale_table = "intermediates/normalization_factors.tsv",  # Produced by an upstream rule in THIS Snakefile.
    output:
        bg_plus_scaled  = "tracks/bedgraphs/{sample}.plus.scaled.bg",
        bg_minus_scaled = "tracks/bedgraphs/{sample}.minus.scaled.bg",
        bw_plus_scaled  = "tracks/bigwigs/{sample}.plus.scaled.bw",
        bw_minus_scaled = "tracks/bigwigs/{sample}.minus.scaled.bw",
    params:
        div_factor = get_div_factor_from_table
    resources:
        mem_mb= 16000,
        runtime= 60
    conda:
        "envs/bedbam_tools.yaml"
    log:
        stdout = "stdout/make_scaled_bigwig/{sample}.out",
        stderr = "stderr/make_scaled_bigwig/{sample}.err"
    threads: 1
    shell:
        """
        # Scale bedGraphs: scaled = unscaled / factor.
        awk -v f="{params.div_factor}" 'BEGIN{{OFS="\t"}} {{ $4 = $4 / f; print }}' {input.bg_plus} \
          > {output.bg_plus_scaled} 2>> {log.stderr}

        awk -v f="{params.div_factor}" 'BEGIN{{OFS="\t"}} {{ $4 = $4 / f; print }}' {input.bg_minus} \
          > {output.bg_minus_scaled} 2>> {log.stderr}

        # Convert scaled bedGraphs -> scaled bigWigs.
        bedGraphToBigWig {output.bg_plus_scaled}  {input.chrom_sizes} {output.bw_plus_scaled}  >> {log.stdout} 2>> {log.stderr}
        bedGraphToBigWig {output.bg_minus_scaled} {input.chrom_sizes} {output.bw_minus_scaled} >> {log.stdout} 2>> {log.stderr}
        """