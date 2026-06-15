# Load necessary packages. 
suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(glue)
})

# Load arguments. 
args = commandArgs(trailingOnly=TRUE)
window_counts = read.csv(args[1], sep = "\t")
sf_path = args[2]
norm_method = args[3]
pval_threshold = as.numeric(args[4])
log2fc_threshold = as.numeric(args[5])
substrate_control = args[6]
tool_dir = args[7]

# Load functions from utils. 
source(glue("{tool_dir}/utils.R"))

# Load the condition names previously defined. 
load("intermediates/conditions.RData")

# Define the columns to keep. 
ctrl_cols  = get_reps(window_counts, ctrl)
trt_cols  = get_reps(window_counts, trt)   
select_cols = c(ctrl_cols, trt_cols)

if (norm_method == "MOR_hier") {
    # Load size factor dataframe and convert to named vector.
    df = read.csv(sf_path, sep = "\t")
    sf = setNames(df$SizeFactor, df$Sample)
} else {
    # Load normalization factor table.
    sf = read.csv(sf_path, sep = "\t")
}

# Run the differential analysis. 
result = process_differential(window_counts, sf, ctrl, trt, select_cols, pval_threshold, log2fc_threshold)
groups=c("IP","Trt") # for PCA analysis

# Save results. 
write.table(result$result_check, file = "tables/all_windows.tsv", sep="\t", row.names = FALSE, quote = FALSE)
write.table(result$sigs, file = "tables/differential_windows.tsv", sep="\t", row.names = FALSE, quote = FALSE)

# For homer analysis. 
write.table(subset(result$sigs, log2FoldChange > 0), file = "intermediates/upreg_windows.tsv", sep="\t", row.names = FALSE, quote = FALSE)
write.table(subset(result$sigs, log2FoldChange < 0), file = "intermediates/downreg_windows.tsv", sep="\t", row.names = FALSE, quote = FALSE)

# Create plots. 
Volcano = volcano_plot(result$result_check, pval_threshold, log2fc_threshold)
FeatureBar = feature_type_bar_plot(result$result_check, pval_threshold, log2fc_threshold)
for_pca = rlog(result$dds)
PCA = plotPCA(for_pca, intgroup=groups)

# Save plots. 
ggsave("plots/Volcano.png", plot = Volcano, width = 6, height = 4, device = "png")
ggsave("plots/FeatureBar.svg", plot = FeatureBar, width = 6, height = 4, device = "svg")
ggsave("plots/PCA.svg", plot = PCA, width = 4, height = 4, device = "svg")

# Create and save deseq plots.
png("plots/MAplot.png", width = 1200, height = 1000, res = 300)
plotMA(result$dds) 
dev.off()
png("plots/DispEst.png", width = 1200, height = 1000, res = 300)
plotDispEsts(result$dds) 
dev.off()

# Create and save results with normalized pseudocounts.
if (norm_method == "MOR_hier") {
    df1 = result$result_check 

    # Find common columns (excluding gene_id, name, and base feature id).
    common_cols = setdiff(intersect(names(df1), names(sf)), c("gene_id", "name", "base_feature_id"))
    
    # Divide each column in df1 by the corresponding value in sf.
    divided = sweep(df1[common_cols], 2, sf[common_cols], `/`)
    divided = ceiling(divided)
    
    # Rename divided columns with "_norm".
    colnames(divided) = paste0(common_cols, "_norm")
    
    # Recombine, gene_id/base_feature_id, unique non-count columns, and normalized values.
    output = cbind(df1[c("gene_id", "base_feature_id")], df1[setdiff(names(df1), c(common_cols, "gene_id", "base_feature_id"))],divided)

    # Save results with the normalized counts. 
    write.table(output, file = "tables/all_windows_norm_counts.tsv", sep="\t", row.names = FALSE, quote = FALSE)
    output_sig = dplyr::filter(output, padj < pval_threshold, abs(log2FoldChange) > log2fc_threshold)
    write.table(output_sig, file = "tables/differential_windows_norm_counts.tsv", sep="\t", row.names = FALSE, quote = FALSE)
} else {
    df1 = result$result_check
    df2 = sf

    # Choose substrate-control ID column.
    if (substrate_control == "INg") {
        id_col = "gene_id"
    } else if (substrate_control == "INf") {
        id_col = "base_feature_id"
    }
    
    # Ensure both have the selected ID column.
    stopifnot(id_col %in% names(df1), id_col %in% names(df2))

    # Reorder df2 to match df1's selected ID.
    df2 = df2[match(df1[[id_col]], df2[[id_col]]), ]

    # Make sure every df1 row found a matching df2 row.
    stopifnot(!any(is.na(df2[[id_col]])))

    # Find shared numeric/count columns, excluding metadata IDs.
    common_cols = setdiff(intersect(names(df1), names(df2)),c("gene_id", "base_feature_id", "name"))

    # Divide element-wise and round up.
    divided = ceiling(df1[common_cols] / df2[common_cols])

    # Rename divided columns with "_norm".
    colnames(divided) = paste0(common_cols, "_norm")

    # Recombine selected ID + unique non-count columns from df1 + normalized counts.
    output = cbind(df1[id_col], df1[setdiff(names(df1), c(common_cols, id_col))], divided)

    # Save results with the normalized counts.
    write.table(output, file = "tables/all_windows_norm_counts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

    # Save significant results with the normalized counts. 
    output_sig = dplyr::filter(output, padj < pval_threshold, abs(log2FoldChange) > log2fc_threshold)
    write.table(output_sig, file = "tables/differential_windows_norm_counts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

# Identify the relevant columns.
ip_trt_cols = grep(paste0("^", trt, "_IP_\\d+_norm$"),names(output_sig),value = TRUE)
ip_ctrl_cols = grep(paste0("^", ctrl, "_IP_\\d+_norm$"),names(output_sig),value = TRUE)
in_trt_cols = grep(paste0("^", trt, "_IN_\\d+_norm$"),names(output_sig),value = TRUE)
in_ctrl_cols = grep(paste0("^", ctrl, "_IN_\\d+_norm$"),names(output_sig),value = TRUE)

# Helper for summing unique features. 
sum_unique_feature_in = function(x, feature_id) {
    keep = !duplicated(feature_id)
    sum(x[keep], na.rm = TRUE)
}

# Sum IPs by gene, sum INs based on substrate control type. 
gene_ip_in = output_sig %>%
    group_by(gene_id, gene_name) %>%
    summarise(
        # Sum IP counts across windows per replicate.
        across(all_of(ip_trt_cols),  ~ sum(.x, na.rm = TRUE)),
        across(all_of(ip_ctrl_cols), ~ sum(.x, na.rm = TRUE)),

        # Aggregate substrate-control counts.
        {
            if (substrate_control == "INg") {
                across(all_of(c(in_trt_cols, in_ctrl_cols)), ~ dplyr::first(.x))
            } else if (substrate_control == "INf") {
                across(
                    all_of(c(in_trt_cols, in_ctrl_cols)),
                    ~ sum_unique_feature_in(.x, .data[["base_feature_id"]])
                )
            } else {
                stop("Unknown substrate_control: ", substrate_control)
            }
        },

        # Calculate the combined fisher P-value. 
        fisher_pval = combine_p_fisher(padj),
        min_padj = {
            m = suppressWarnings(min(padj, na.rm = TRUE))
            if (is.infinite(m)) NA_real_ else m
        },

        # Calculate the log fold-change at whichever window has the minimum p-value. 
        lfc_at_min_p = {
            p = padj
            i = which.min(replace(p, is.na(p), Inf))
            if (all(is.na(p))) NA_real_ else log2FoldChange[i]
        },

        n_sites = n(),
        .groups = "drop"
    )

# Average IP and IN across replicates within each condition.
gene_ip_in = gene_ip_in %>%
    mutate(IP_trt_mean = if (length(ip_trt_cols)  > 0) rowMeans(across(all_of(ip_trt_cols)),  na.rm = TRUE) else NA_real_,
            IP_ctrl_mean = if (length(ip_ctrl_cols) > 0) rowMeans(across(all_of(ip_ctrl_cols)), na.rm = TRUE) else NA_real_,
            IN_trt_mean = if (length(in_trt_cols)  > 0) rowMeans(across(all_of(in_trt_cols)),  na.rm = TRUE) else NA_real_,
            IN_ctrl_mean = if (length(in_ctrl_cols) > 0) rowMeans(across(all_of(in_ctrl_cols)), na.rm = TRUE) else NA_real_)

# Compute new, total log-fold change.
pseudo = 0.2 

# Calculate total log fold-change (sort of like the average log fold-change across all sites). 
gene_ip_in = gene_ip_in %>%
    mutate(ratio_trt  = (IP_trt_mean  + pseudo) / (IN_trt_mean  + pseudo), ratio_ctrl = (IP_ctrl_mean + pseudo) / (IN_ctrl_mean + pseudo),
           total_log2FC = log2(ratio_trt / ratio_ctrl))

# Subset to just the columns of interest. 
gene_summary = select(gene_ip_in, c('gene_id', 'gene_name', 'total_log2FC', 'fisher_pval','min_padj','lfc_at_min_p','n_sites')) 

# Save the gene level data. 
write.table(gene_summary, file = "tables/differential_genes.tsv", sep="\t", row.names = FALSE)

# Min-padj approach.
dot_plot_1 = plot_top_genes(gene_summary, "min_padj", "lfc_at_min_p", score_label = "Minimum log10(p-value)")
ggsave("plots/genes_by_top_site.svg", dot_plot_1, width = 8, height = 9, device = "svg")

# Fisher p-value approach.
dot_plot_2 = plot_top_genes(gene_summary, "fisher_pval", "total_log2FC", score_label = expression(-log[10](Fisher~p)))
ggsave("plots/genes_by_all_sites.svg", dot_plot_2, width = 8, height = 9, device = "svg")

# Make a color-coded BED.

# Define track name and description. 
track_name = "Flipper Differential Windows"
track_desc = "Positive=red, Negative=blue"

# Coerce coords.
result$sigs[["start"]] = as.integer(result$sigs[["start"]])
result$sigs[["end"]] = as.integer(result$sigs[["end"]])

# Coerce LFC to numeric.
result$sigs[["log2FoldChange"]] = suppressWarnings(as.numeric(result$sigs[["log2FoldChange"]]))

# Map the abolute log foldchange to the BED score 
score_fun = function(lfc) {
    s = round(1000 * (1 - exp(-abs(lfc))))
    s[is.na(s)] = 0
    pmax(0, pmin(1000, as.integer(s)))
}

# Assign colours assignment.
rgb = rep("128,128,128", nrow(result$sigs))  # default gray
rgb[!is.na(result$sigs[["log2FoldChange"]]) & result$sigs[["log2FoldChange"]] > 0] = "255,0,0"
rgb[!is.na(result$sigs[["log2FoldChange"]]) & result$sigs[["log2FoldChange"]] < 0] = "0,0,255"

# Set strand to blank (.)
strand = rep(".", nrow(result$sigs))

# Thick start/end (looks better).
thickStart = result$sigs[["start"]]
thickEnd = result$sigs[["end"]]

# Create the bed. 
bed = data.frame(chrom = result$sigs[["chr"]], chromStart = result$sigs[["start"]], chromEnd = result$sigs[["end"]],
                  name = as.character(seq_len(nrow(result$sigs))), score = score_fun(result$sigs[["log2FoldChange"]]),
                  strand = strand, thickStart = thickStart, thickEnd = thickEnd, itemRgb = rgb,
                  stringsAsFactors = FALSE)

# Save the bed to file. 
con = file("tables/differential_windows.bed", open = "wt")
writeLines(sprintf('track name="%s" description="%s" itemRgb="On" visibility=3', track_name, track_desc), con = con)
write.table(bed, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)

