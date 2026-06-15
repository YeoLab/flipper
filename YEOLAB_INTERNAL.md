# NOTE TOO SELF: I should really save this for the time when I actually have the module set up so I do not need to walk people through all of this nonsense. 

# YEO-LAB internal example 
Hello. This is a short example for running Flipper as a member of the Yeo-lab partition on TSCC.

## Load up an interactive node.

```
srun -N 1 -c 1 -t 4:00:00 -p gold -q hcg-csd792 -A csd792 --mem 4G --pty /bin/bash
```

## Create a folder in scratch to save the output.
Run the command below (Replacing YOUR_USERNAME with your TSCC username) to create a folder to save the output of the 

```
mkdir /tscc/lustre/ddn/scratch/YOUR_USERNAME/Flipper_test
cd /tscc/lustre/ddn/scratch/YOUR_USERNAME/Flipper_test
```

## Load up the skipper module.

```
module load skipper
```

Yes, you load the Skipper module to run Flipper, as they have identical dependencies. In the future, we will most likely 

**IMPORTANT**: The only modules that should be loaded are the default modules (the ones that automatically load whenever you go on TSCC) and Skipper (Loading Skipper also automatically loads the singularity module). No other modules should be loaded, as this may confuse snakemake as to which environment to run the code in, causing an error. You can check what modules you have loaded using ```module list ```

## Preparing a new Flipper run:

### Config file

Flipper includes a range of parameters to accommodate user preferences and diverse biological questions. These settings are defined in a config.yaml file, which is supplied to Snakemake during execution. An example configuration file (Flipper.yaml) is provided in this repository, and detailed descriptions of all available options are listed below.

| Command    | Description |
| -------- | ------- |
| CONTROL  | The name used in the "Sample" column of your skipper manifest that corresponds to whatever you consider to be the control   |
| TREATMENT  | The name used in the "Sample" column of your skipper manifest that corresponds to whatever you consider to be the treatment   |
| SKIPPER_CONFIG | The path to the config file used by skipper |
| SKIPPER_OUTPUT   | The path to the folder where skipper saved its output    |
| OUTPUT | The path to the folder where the Flipper results will be saved |
| LOG2FC_THRESHOLD | Absolute log2FC threshold for filtering |
| PVAL_THRESHOLD   | P-value threshold for filtering    |
| NORMALIZATION | Currently supported normalization techniques include hierarchical median of ratios ("MOR_hier") and hierarchical EDAseq (EDA_hier) |
| EDA_GC| The GC normalization method used by EDAseq (See [EDAseq manual](https://bioconductor.org/packages/devel/bioc/manuals/EDASeq/man/EDASeq.pdf)). NOTE: Ignored if "EDA_hier" was not chosen for NORMALIZATION |
| EDA_LENGTH| The length normalization method used by EDAseq (See [EDAseq manual](https://bioconductor.org/packages/devel/bioc/manuals/EDASeq/man/EDASeq.pdf])). NOTE: Ignored if "EDA_hier" was not chosen for NORMALIZATION |
| EDA_SAMPLE| The Sample normalization method used by EDAseq (See [EDAseq manual](https://bioconductor.org/packages/devel/bioc/manuals/EDASeq/man/EDASeq.pdf])). NOTE: Ignored if "EDA_hier" was not chosen for NORMALIZATION |
| GO | Optional: Run a GO analysis on genes found to have significant up/down regulation (true/false) |
| SPECIES | Optional: Used for selecting accurate terms during GO analysis. species must following the naming convention used by the annotation hub package in R (e.g. "Homo sapiens", "Mus musculus", "Danio rerio")
| HOMER | Optional: Generate homer motif output for differential windows (true/false) |
| METADENSITY | Optional: Generate protein binding densities for differential windows (true/false) |
| RNASEQ | Optional: Path to RNASEQ DESeq2 result for a combined gene expression differential window analysis (file path).|

NOTE ON NORMALIZATION: Normalizing eCLIP data across seperate treatments is complicated, and defining a singular best normalization method is difficult. Generally speaking, MOR_hier normalization is sufficient for most use cases, while EDAseq_hier is useful if you want extra control over confounding factors like GC content and length. Please see the "Choosing a Normalization Method" section for more details.

### Choosing a normalization method: 

- Median of Ratios hierarchical (MOR_hier) is a simple and robust normalization option, making it ideal for most differential RBP binding analyses. However, it does have some caveats, such as not correcting for GC-bias or length (please note that identical regions are compared across conditions, making the effect of these biases less significant than if one was comparisoning binding between different regions). 
- EDAseq hierarchical (EDA_hier) normalization utilizes more advanced normalization techniques that attempt to control for additional factors such as GC and length bias. As such, EDAseq offers several normalization techniques for each of its three normalization steps, GC normalization, length normalization, and traditional cross sample normalization. For more information on EDAseq normalization, please see the [original publication](https://pubmed.ncbi.nlm.nih.gov/22177264/) and the [vignette](https://bioconductor.org/packages/release/bioc/vignettes/EDASeq/inst/doc/EDASeq.html). We have personally found that using the full quantile normalization (full) for all options works best.
    - Note that we cannot guaruntee that every combination of EDAseq normalization options will work for all datasets.

For more information on normalization techniques, including Flippers hierarchical IP normalization method, Please see the [preprint](https://doi.org/10.64898/2026.03.13.711628)

## Flipper Output

Flipper generates a variety of output files that are organized into different folders.

### Tables

Flipper creates several tab-delimited tables that summarize the results of the differential binding analysis:
- differential_windows.tsv
    - Contains all windows found to exhibit significantly differential binding following treatment (adjusted p-value ≤ 0.05). The table is sorted by adjusted p-value and includes useful metadata for each significant site.
- all_windows.tsv
    - Contains the same columns as differential_windows.tsv, but includes every tested window, not just the significant ones.
- differential_windows_norm_counts.tsv and all_windows_norm_counts.tsv
    - Identical to the corresponding tables above, but with raw read counts replaced by normalized counts (i.e., counts divided by size factors). This allows users to better visualize why particular windows were or were not classified as significant.
    - differential_genes.tsv
    - Summarizes significant results at the gene level. For each gene, it reports:
        - The number of differential windows contained within the gene.
        - The mean log fold change across all differential windows.
        - Fisher’s combined p-value (aggregate significance for the gene).
        - The minimum p-value observed among the gene’s windows.
        - The log fold change corresponding to the window with the minimum p-value.

### Plots

Flipper produces a variety of plots to help visualize and interpret differential binding results:
- Volcano plots — Show the distribution of log fold changes vs. p-values for all tested windows.
- Feature type barplots — Show how many windows were found to be upregulated or downregulated within each feature type (e.g., introns, CDS regions).
- Top-20 gene dotplots — Summarize read-level signal for the top 20 genes with the strongest differential binding.
- PCA plots (from DESeq2) — Useful for checking that within-group replicates cluster more closely than between-group replicates.
- MA plots (from DESeq2) — Visualize the dispersion and effect sizes of windows.
- Dispersion estimates (from DESeq2) — Show the dispersion–mean relationship estimated by DESeq2.
- Gene Ontology enrichment — Highlights enriched GO terms for genes with upregulated or downregulated windows.
- Metadensity - Shows the relative 

For more information on metadensity, please see the [original publication](https://pubmed.ncbi.nlm.nih.gov/36388152/)

### HOMER Motif Analysis

Flipper automatically runs HOMER motif analysis on significantly upregulated and downregulated windows. This helps identify whether specific sequence motifs are over-represented within affected regions. For more details on these outputs, please see [HOMER](http://homer.ucsd.edu/homer/motif/)

### Intermediates

Flipper saves several intermediate datasets produced during differential site calling. These files can be useful for custom analyses, but most users can ignore them for routine use.

### Tracks

Skipper outputs bedgraphs and bigwigs for each replicate. These outputs have been scaled using the same normalization factors used by Flipper during analysis, and can thus be helpful to understand how Flipper identifies differential regions. 

For unscaled tracks, please see the Skipper output. 

# Troubleshooting

1. Log files.
Flipper generates 3 types of log files. The first 2 types can be found within the `stderr/` and `stdout/` folders of your WORKDIR. These log files generally should contain all the information needed for debugging. 

However, in some cases additional information from snakemake may be necessary, in which cases users are encouraged to investigate the log files in `WORKDIR/.snakemake/slurm_logs`. These log files are organized by rules and contain additional information on the snakemake run.

2. TODO: troubleshooting section for if the initial conda installation stalls. 

3. Problems with EDAseq_hier normalization.
If you are using EDAseq_hier normalization and you find your snakemake runs are reguarly breaking at the normalization rule, it is likely some problem with the combination of EDAseq normalization inputs you are giving to Flipper. Our first reccomendation is to try the run with MOR_hier normalization to confirm that the problem is with EDAseq. We then reccomend trying out a few different combinations of EDAseq normalization parameters to see if you can find one that works. Again, as we cannot guaruntee that EDAseq is suitable for all datasets, defaulting to MOR_hier my be necessary. 

4. Jobs dying with no explanation.
If you observe that many of your jobs are dying without any explanation (e.g. mostly blank files in WORKDIR/stderr, unhelpful error messages in WORKDIR/.snakemake/slurm_logs such as "Killed"), and these jobs are occuring on the same node according to WORKDIR/stdout, then it is likely that this is the result of problematic nodes on your cluster. I would reccomend taking whichever nodes were used for the failed jobs and excluding them from the analysis by adding the following lines to the slurm extra command within your profile like so:

  ```yaml
  slurm_extra: "--exclude=YOUR_NODE"
  ```

This is also the common cause of many timeout errors, as Skipper generally provides significantly more than enough time for all rules. 

5. Monitoring pipeline. 
`squeue -u $USER -o "%.18i %.10P %.20j %.10u %.2t %.10M %.6D %.20R %.80k"` will show currently active jobs (helpful to see if certain rules are getting stuck.)