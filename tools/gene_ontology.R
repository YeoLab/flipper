################### Initial Setup ##################

# Load necessary packages. 
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(AnnotationHub)
    library(AnnotationDbi)
    library(tidyverse)
    library(cowplot)
    library(ggplot2)
    library(enrichplot)
    library(glue)
})

# Load up the arguments passed from snakemake.
args = commandArgs(trailingOnly=TRUE)

# Divide arguments into proper names. 
species = args[1] 
diff_result = args[2]
gene_wise_IN = args[3]
output = args[4]
tool_dir = args[5]

# Load up functions from utils. 
source(glue("{tool_dir}/utils.R"))

################### Setup annotation database ########################

#Load up annotation hub
hub = AnnotationHub()
orgdbs = query(hub, c("OrgDb", species))
orgdb = orgdbs[[length(orgdbs)]]

# Extract just the relationship between the go term and the ensembl ID. 
go_annot = AnnotationDbi::select(orgdb, keys = keys(orgdb, keytype = "ENSEMBL"), columns = c("GO"), keytype = "ENSEMBL")
term2gene = go_annot[, c("GO", "ENSEMBL")]

################ Load the data ###################

# Load up the differential genes. 
diff_genes = read.csv(diff_result, sep = "\t")

# Split differential genes into upregulated and downregulated genes. 
diff_genes_up = subset(diff_genes, total_log2FC > 0)$gene_id
diff_genes_down = subset(diff_genes, total_log2FC < 0)$gene_id

# Load up the gene_wise_IN values
gene_wise_IN = read.csv(gene_wise_IN, sep = "\t")

# Filter out all genes with extremely low expression, as they should not contribute to the background. 
in_ip_cols = grep("_IN_|_IP_", names(gene_wise_IN), value = TRUE)
num_cols = length(in_ip_cols)
INg_filtered = gene_wise_IN %>% mutate(total_expr = rowSums(pick(all_of(in_ip_cols)), na.rm = TRUE)) %>% filter(total_expr > num_cols)

# Extract the column out into a vector.
gene_vec = INg_filtered$gene_id

# Split multi-gene entries into separate entries.
gene_vec_split = str_split(gene_vec, ":", simplify = FALSE) %>% unlist()

# Remove everything after "."
gene_vec_clean = str_remove(gene_vec_split, "\\..*")

# Ensure uniqueness.
gene_vec_unique = unique(gene_vec_clean)

############### Run the GO analysis #############

# Make sure to strip version suffixes from both sets if needed.
gene_vec_unique = sub("\\..*", "", gene_vec_unique)
diff_genes_up  = sub("\\..*", "", diff_genes_up)
diff_genes_down  = sub("\\..*", "", diff_genes_down)

# Run enrichment analysis for upregulated. 
compMF_up = enrichGO(gene = diff_genes_up, universe = gene_vec_unique, OrgDb = orgdb, keyType = "ENSEMBL",
                      ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

# Run enrichment analysis for downregulated. 
compMF_down = enrichGO(gene = diff_genes_down, universe = gene_vec_unique, OrgDb = orgdb, keyType = "ENSEMBL",         
                        ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

################ Plotting #################

# Ensures all y-axis labels are readable. 
wrap_y = function(p, width = 40, left_margin = 140) {
    p + scale_y_discrete(labels = function(x) str_wrap(x, width)) +
        coord_cartesian(clip = "off") + theme(plot.margin = margin(10, 20, 10, left_margin))
}

# Create placeholder if there are no enriched terms.              
empty_panel = function(title = "No enriched terms") {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = title, size = 5, fontface = "bold") +
        theme_void() + theme(plot.margin = margin(10, 20, 10, 20))
}

# Build dotplot.                        
make_dotplot = function(x, title, show = 10) {
    has_rows = tryCatch({
        df = as.data.frame(x)
        is.data.frame(df) && nrow(df) > 0
    }, error = function(...) FALSE)

    if (!has_rows) return(empty_panel(title)) # Checks to make sure it is not empty. 

    p = dotplot(x, showCategory = show) + ggtitle(title) +
        theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              axis.text  = element_text(size = 13), axis.title = element_text(size = 13))

    wrap_y(p, width = 38, left_margin = 160)
}

# Build panels individually, ensuring no text gets cut off. 
p_up = make_dotplot(compMF_up, "GO Enrichment: Upregulated Genes", show = 10)
p_down = make_dotplot(compMF_down, "GO Enrichment: Downregulated Genes", show = 10)

# Combine panels with a bit of buffer space. 
combined = plot_grid(p_up, NULL, p_down, ncol = 3, rel_widths = c(1, 0.01, 1), align = "h")

# Save the combined plot.
ggsave(filename = output, plot = combined, device = "svg", width = 20, height = 6.5)