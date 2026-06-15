# Flipper
![Flipper cartoon](documents/logo.png)

For finding when binding flips.

A DESeq2-based method for identifying significant differential binding sites using eCLIP data

See preprint at: [https://doi.org/10.64898/2026.03.13.711628](https://doi.org/10.64898/2026.03.13.711628)

# Set up
## Installing and Running Skipper: 

Flipper acts as an extension to the Skipper program previously developed by the Yeo Lab. As such, Skipper output is required for running Flipper. Skipper installation instructions can be found [here](https://github.com/YeoLab/skipper). 

Please note that Flipper and Skipper operate under similar frameworks. As such, many of Flipper's installation instructions are identical/redundant with the installation instructions for Skipper. That being said, all Flipper installation instructions are still included for completeness. 

## Installing Flipper:

1. **Clone the repository**  
   ```bash
   git clone https://github.com/YeoLab/flipper.git
   cd flipper
   ```

2. **Install Conda (if not already installed)**  
   The example below shows how to install Miniconda on Linux (64-bit). For detailed instructions on other systems, see the [official installation guide](https://www.anaconda.com/docs/getting-started/miniconda/install).  

   ```bash
   curl -L -O "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

3. **Create a Snakemake environment**  
   Once Conda is installed, create an environment with Snakemake version `9.12.0` to run Flipper:  

   ```bash
   conda create -n snakemake9 snakemake=9.12.0
   ```

All other required packages and environments will be installed automatically by Snakemake and Conda the first time you run Flipper.

## Configuring Your Snakemake Profile

Snakemake profiles allow you to supply additional arguments without cluttering the command line.  
An example profile is provided at:`profiles/example_basic/config.yaml`

This profile is configured for running Flipper on a single-node machine. The only required change is to specify a path for saving Conda environments (choose any location on your machine with sufficient storage space).

## Running Flipper on HPCs

Flipper's core functionality is lightweight enough to run on most personal computers. However, as an extension of the much heavier Skipper program, Flipper is also optimized to run on the same high performance clusters (HPCs) usually needed for Skipper. 

### Cluster Executor Setup

To run Flipper on HPCs, you must install a cluster executor plugin.  
The example below demonstrates installation of the **SLURM** executor (a widely used workload manager).  
Other executor options are listed in the [Snakemake plugin catalog](https://snakemake.github.io/snakemake-plugin-catalog/index.html).

```bash
conda activate snakemake9
conda install snakemake-executor-plugin-slurm=1.4.0
```

### Adjusting Your Profile

After installing the executor plugin, you must adjust your Snakemake profile.  
An example profile is provided in:

```
profiles/example_slurm/config.yaml
```

- **CONDA prefix** Do not forget to change this to a path on your cluster. 
- **Slurm Account and partition:** You must enter your own account and partition information.
- **Cluster-specific options:** Some systems require additional details. For example:  

  ```yaml
  slurm_extra: "--qos=YOUR_QOS"
  ```

  In general, if something is required in your `srun` or `sbatch` commands, it may also need to be added to `slurm_extra`.


## Minimal Example. 

This section provides a small example run of Flipper. The example uses a heavily modified version of Skipper output that includes only the inputs required for Flipper and has been subset to chromosome 22. The underlying Skipper run is derived from the NONO dataset described in [https://doi.org/10.64898/2026.03.13.711628](https://doi.org/10.64898/2026.03.13.711628), which was originally generated in [https://doi.org/10.1038/s41589-023-01270-0](https://doi.org/10.1038/s41589-023-01270-0).

This example assumes that you are working on a Linux-based system with Slurm set up and have already gone through all installation steps above (including adjusting the example profile).

1. **Edit config files**
   Open both the `example/example_flipper_config.yaml` and the `example/example_skipper_config.yaml` config files using any text editor and change every instance of `/path/to/your/Flipper` to the **absolute** path to the Flipper directory you just cloned. Also, change every instance of `/path/to/save/output` with the **absolute** path to whatever location you want to save your Flipper outputs to.

Please note that this example has both the "MAKE_BIGWIGS" and "HOMER" options set to false. This is done to minimize processing time and input size for this example. For more information on these options, please see the [Config file](#config-file) section. 

2. **Run Flipper**  
   ```bash
   unset SLURM_JOB_ID # required if running on an interactive node, which is recommended
   
   snakemake -ks Flipper.py --configfile example/example_flipper_config.yaml --profile ../profiles/tscc2_snakemake9
   ```

NOTE: The first run of Flipper needs to set up all of the necessary conda environments via snakeconda. As such, this initial Flipper run will be quite slow, but subsequent runs will be much faster.

NOTE: If difficulties arise while running this example (or any run of Flipper) please see the [Troubleshooting](#Troubleshooting) section and/or open an issue. 

## Preparing a new Flipper run:

### Config file

Flipper includes a range of parameters to accommodate user preferences and diverse biological questions. These settings are defined in a config.yaml file, which is supplied to Snakemake during execution. An example configuration file (`example/example_flipper_config.yaml`) is provided in this repository, and detailed descriptions of all available options are listed below.

| Command    | Description |
| -------- | ------- |
| SKIPPER_CONFIG | The path to the config file used by Skipper |
| SKIPPER_OUTPUT   | The path to the directory where Skipper saved its output    |
| WORKDIR | The path to the directory where the Flipper results will be saved |
| TOOLS | The path to the `tools` directory from this repository |
| CONTROL  | The name used in the "Sample" column of your Skipper manifest that corresponds to whatever you consider to be the control   |
| TREATMENT  | The name used in the "Sample" column of your Skipper manifest that corresponds to whatever you consider to be the treatment   |
| LOG2FC_THRESHOLD | Absolute log2FC threshold for filtering |
| PVAL_THRESHOLD   | P-value threshold for filtering |
| SUBSTRATE_CONTROL | Specify whether to sum up IN values across genes ("INg") or genomic features ("INf") |
| NORMALIZATION | Currently supported normalization techniques include hierarchical median of ratios ("MOR_hier") and hierarchical EDASeq (EDA_hier) |
| SPECIES | Used for selecting accurate terms during GO analysis. Species must follow the naming convention used by the annotation hub package in R (e.g. "Homo sapiens", "Mus musculus", "Danio rerio") |
| GC_METHOD | The GC normalization method used by EDASeq (See [EDASeq manual](https://bioconductor.org/packages/devel/bioc/manuals/EDASeq/man/EDASeq.pdf)). NOTE: Ignored if "EDA_hier" was not chosen for NORMALIZATION |
| LENGTH_METHOD | The length normalization method used by EDASeq (See [EDASeq manual](https://bioconductor.org/packages/devel/bioc/manuals/EDASeq/man/EDASeq.pdf])). NOTE: Ignored if "EDA_hier" was not chosen for NORMALIZATION |
| BETWEEN_METHOD | The sample normalization method used by EDASeq (See [EDASeq manual](https://bioconductor.org/packages/devel/bioc/manuals/EDASeq/man/EDASeq.pdf])). NOTE: Ignored if "EDA_hier" was not chosen for NORMALIZATION |
| MAKE_BIGWIGS | Optional: Generate normalized bigWig files using the same normalization scheme selected for Flipper. For unnormalized bigWig files, please use the Skipper output |
| HOMER | Optional: Generate HOMER motif output for differential windows (True/False) |
| METADENSITY | Optional: Generate protein binding densities for differential windows (True/False) |
| FEATURES | Optional: A list of feature types to use for the metadensity plots. Required if METADENSITY is True. Options include: 'five_prime_UTR', 'first_exon', 'exon', 'intron', 'last_exon', 'CDS', and 'three_prime_UTR'|

NOTE ON NORMALIZATION: Normalizing eCLIP data across seperate treatments is complicated, and defining a singular best normalization method is difficult. Generally speaking, MOR_hier normalization is sufficient for most use cases, while EDASeq_hier is useful if you want extra control over confounding factors like GC content and length. Please see the "Choosing a Normalization Method" section for more details.

### Choosing a normalization method: 

- Median of Ratios hierarchical (MOR_hier) is a simple and robust normalization option, making it ideal for most differential RBP binding analyses. However, it does have some caveats, such as not correcting for GC-bias or length (please note that identical regions are compared across conditions, making the effect of these biases less significant than if one was comparing binding between different regions). 
- EDASeq hierarchical (EDA_hier) normalization utilizes more advanced normalization techniques that attempt to control for additional factors such as GC and length bias. As such, EDASeq offers several normalization techniques for each of its three normalization steps, GC normalization, length normalization, and traditional cross-sample normalization. For more information on EDASeq normalization, please see the [original publication](https://pubmed.ncbi.nlm.nih.gov/22177264/) and the [vignette](https://bioconductor.org/packages/release/bioc/vignettes/EDASeq/inst/doc/EDASeq.html). We have personally found that using the full quantile normalization (full) for all options works best.
    - Note that we cannot guarantee that every combination of EDASeq normalization options will work for all datasets.

For more information on normalization techniques, including Flipper's hierarchical IP normalization method, please see the [preprint](https://doi.org/10.64898/2026.03.13.711628)

## Flipper Output

Flipper generates a variety of output files that are organized into different folders.

### Tables

Flipper creates several tab-delimited tables that summarize the results of the differential binding analysis:
- differential_windows.tsv
    - Contains all windows found to exhibit significantly differential binding following treatment (adjusted p-value ≤ 0.05). The table is sorted by adjusted p-value and includes useful metadata for each significant site.
- all_windows.tsv
    - Contains the same columns as differential_windows.tsv, but includes every tested window, not just the significant ones.
- differential_windows_norm_counts.tsv and all_windows_norm_counts.tsv
    - Identical to the corresponding tables above, but with raw read counts replaced by normalized counts (i.e., counts adjusted using the selected normalization method). This allows users to better visualize why particular windows were or were not classified as significant.
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
- Metadensity — Shows the relative information of the CLIP data across user specified regions. For more information on metadensity, please see the [original publication](https://pubmed.ncbi.nlm.nih.gov/36388152/)

### HOMER Motif Analysis

Flipper automatically runs HOMER motif analysis on significantly upregulated and downregulated windows. This helps identify whether specific sequence motifs are over-represented within affected regions. For more details on these outputs, please see [HOMER](http://homer.ucsd.edu/homer/motif/).

### Intermediates

Flipper saves several intermediate datasets produced during differential site calling. These files can be useful for custom analyses, but most users can ignore them for routine use.

### Tracks

Skipper outputs bedGraphs and bigWigs for each replicate. These outputs have been scaled using the same normalization factors used by Flipper during analysis, and can thus be helpful to understand how Flipper identifies differential regions. 

For unscaled tracks, please see the Skipper output. 

# manuscript figures

This folder contains notebooks specifically designed for creating the figures for the Flipper publication. These files are included to ensure total reproducibility of the results included in the original Flipper manuscript. Most users of Flipper can ignore the contents of this folder. 

# Troubleshooting

1. Log files.
Flipper generates 3 types of log files. The first 2 types can be found within the `stderr/` and `stdout/` folders of your WORKDIR. These log files generally should contain all the information needed for debugging. 

However, in some cases additional information from Snakemake may be necessary, in which cases users are encouraged to investigate the log files in `WORKDIR/.snakemake/slurm_logs`. These log files are organized by rules and contain additional information on the Snakemake run.

2. TODO: troubleshooting section for if the initial conda installation stalls. 

3. Problems with EDASeq_hier normalization.
If you are using EDASeq_hier normalization and you find your Snakemake runs are reguraly breaking at the normalization rule, it is likely that there are some problems with the combination of EDASeq normalization inputs you are giving to Flipper. Our first recommendation is to try the run with MOR_hier normalization to confirm that the problem is with EDASeq. We then recommend trying out a few different combinations of EDASeq normalization parameters to see if you can find one that works. Again, as we cannot guarantee that EDASeq is suitable for all datasets, defaulting to MOR_hier may be necessary. 

4. Jobs dying with no explanation.
If you observe that many of your jobs are dying without any explanation (e.g. mostly blank files in `WORKDIR/stderr`, unhelpful error messages in `WORKDIR/.snakemake/` slurm_logs such as "Killed"), and these jobs are occurring on the same node according to `WORKDIR/stdout`, then it is likely that this is the result of problematic nodes on your cluster. I would recommend taking whichever nodes were used for the failed jobs and excluding them from the analysis by adding the following lines to the slurm extra command within your profile like so:

  ```yaml
  slurm_extra: "--exclude=YOUR_NODE"
  ```
This is also the common cause of many timeout errors, as Flipper generally provides significantly more than enough time for all rules. 

5. Monitoring pipeline. 
`squeue -u $USER -o "%.18i %.10P %.20j %.10u %.2t %.10M %.6D %.20R %.80k"` will show currently active jobs (helpful to see if certain rules are getting stuck.)