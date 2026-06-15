# Load necessary packages. 
suppressPackageStartupMessages({
    library(tidyverse)
    library(DESeq2)
    library(data.table)
    library(EDASeq)
    library(glue)
    library(purrr)
})

# Load up the arguments passed from snakemake. 
args = commandArgs(trailingOnly=TRUE)

# Divide arguments into proper names. 
method = args[1]
sample_ctrl = args[2] 
sample_trt = args[3] 
skipper_dir = args[4]
anno_path = args[5]
tool_dir = args[6]
gc_method = args[7]
length_method = args[8]
between_method = args[9]
substrate_control = args[10]
output = "intermediates/normalization_factors.tsv"

# Load up functions from utils. 
source(glue("{tool_dir}/utils.R"))
load("intermediates/conditions.RData")

# Calculate genewise IN by summing over gene name in the count data. 
count_data = read.csv(paste0("intermediates/all_count_data.tsv"), sep = '\t')
gene_wise = get_substrate_level(count_data, ctrl, trt, "INg")
gene_wise_IN = gene_wise$substrate_wise_IN
gene_wise_IP = gene_wise$substrate_wise_IP
write.table(gene_wise_IN, "intermediates/gene_wise_IN.tsv", sep="\t", row.names = FALSE)
write.table(gene_wise_IP, "intermediates/gene_wise_IP.tsv", sep="\t", row.names = FALSE)

# Load up the condition names previously defined. 
load("intermediates/conditions.RData")

# A weak filter is required to avoid potential errors from many zero count rows. 
ctrl_cols  = get_reps(count_data, ctrl)  
trt_cols  = get_reps(count_data, trt)  
just_counts = c(ctrl_cols, trt_cols)
filtered_counts = count_data[rowMeans(count_data[just_counts]) > 1,]

# Load up the reproducible windows for each treatment group. 
IP_ctrl = get_window_counts(filtered_counts, skipper_dir, sample_ctrl, ctrl)
IP_trt = get_window_counts(filtered_counts, skipper_dir, sample_trt, trt)
window_ctrl = IP_ctrl$window_counts
window_trt = IP_trt$window_counts
background_ctrl = IP_ctrl$background_counts
background_trt = IP_trt$background_counts

# Get the rep counts for later.
reps_trt = as.integer(sub(".*_(IN|IP)_", "", grep(paste0("^", trt, "_"), names(filtered_counts), value = TRUE))) |> unique()
reps_ctrl = as.integer(sub(".*_(IN|IP)_", "", grep(paste0("^", ctrl, "_"), names(filtered_counts), value = TRUE))) |> unique()

if (method == "MOR_hier") {

    # Calculate size factors for the input. 
    IN_sf = get_sf_MOR(gene_wise_IN[, setdiff(names(gene_wise_IN), c("gene_name", "base_feature_id","length", "gc", "name", "gene_id"))])
    
    # Calculate norm factors for the treatment and ctrl IPs seperately. 
    IP_ctrl_sf = get_sf_MOR(window_ctrl[, setdiff(names(window_ctrl), c("gene_name","base_feature_id","length", "gc", "name", "gene_id"))])
    IP_trt_sf = get_sf_MOR(window_trt[, setdiff(names(window_trt), c("gene_name","base_feature_id","length", "gc", "name", "gene_id"))])

    # Combine into one vector.
    temp_norm_factors = c(IP_trt_sf, IP_ctrl_sf, IN_sf)

    # Build the vector of column names for control IP replicates
    ctrl_cols = names(IP_ctrl_sf)
    
    # Same for treatment IP replicates
    trt_cols  = names(IP_trt_sf)
    
    # Keep only background windows present in both conditions. 
    background_both = inner_join(background_ctrl, background_trt, by = "name", suffix = c("_ctrl", "_trt"))
    
    # Sum them as new columns.
    background_both[[paste0(ctrl, "_IP_sum")]] = rowSums(background_both[paste0(ctrl_cols)], na.rm = TRUE)
    background_both[[paste0(trt, "_IP_sum")]]  = rowSums(background_both[paste0(trt_cols)],  na.rm = TRUE)
    
    # Collapse background windows to genes.
    background_gene = background_both %>%
        group_by(gene_id_ctrl) %>%
        summarise(
            across(all_of(c(paste0(ctrl, "_IP_sum"), paste0(trt, "_IP_sum"))), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop"
        ) %>%
        dplyr::rename(gene_name = gene_id_ctrl)
    
    # Size factors for average.
    ctrl_sum_col = paste0(ctrl, "_IP_sum")
    trt_sum_col  = paste0(trt, "_IP_sum")
    IP_backg_sf = get_sf_MOR(dplyr::select(background_gene, gene_name, all_of(c(ctrl_sum_col, trt_sum_col))))

    # Combine the rep specific size factors with the average background size factors. 
    norm_factors = combine_ip_size_factors(temp_norm_factors, IP_backg_sf, ctrl, trt)

    # Save the size factors. 
    write.table(data.frame(Sample = names(norm_factors), SizeFactor = norm_factors),
              output,row.names = FALSE, sep = "\t")
} else if (method == "EDA_hier") { 
    
    # Collect different input values based on the use of substrate control.  
    if (substrate_control == "INf") {
        feature_wise = get_substrate_level(count_data, ctrl, trt, "INf")
        EDA_IN = feature_wise$substrate_wise_IN
        id_col = "base_feature_id"
    } else {
        EDA_IN = gene_wise_IN
        id_col = "gene_id"
    }
    
    # Drop all complete 0 rows, as they cause errors when given to EDAseq. 
    EDA_IN = drop_0_rows(EDA_IN)

    # Drop all rows with zeros in either of these sets (0s are still a problem for EDA, but we want all drops to be identical)
    filtered = drop_0_rows_joint_universe(window_ctrl, background_ctrl, window_trt, background_trt)

    # Re-extract the filtered results. 
    window_ctrl = filtered$window_ctrl
    background_ctrl = filtered$background_ctrl
    window_trt = filtered$window_trt
    background_trt = filtered$background_trt
    
    # Calculate size factors for the input. 
    IN_sf = get_sf_EDA(c(ctrl,trt), EDA_IN, id_col, length(reps_ctrl), gc_method, length_method, between_method)
    
    # Calculate size factors for the IP. 
    IP_signal_sf_ctrl = get_sf_EDA(c(ctrl), dplyr::select(window_ctrl, -all_of(id_col)), "name", length(reps_ctrl),
                                  gc_method, length_method, between_method)
    IP_noise_sf_ctrl = get_sf_EDA(c(ctrl), dplyr::select(background_ctrl, -all_of(id_col)), "name", length(reps_ctrl),
                                 gc_method, length_method, between_method)
    IP_signal_sf_trt = get_sf_EDA(c(trt), dplyr::select(window_trt, -all_of(id_col)), "name",  length(reps_trt),
                                  gc_method, length_method, between_method)
    IP_noise_sf_trt = get_sf_EDA(c(trt), dplyr::select(background_trt, -all_of(id_col)), "name",  length(reps_trt),
                                 gc_method, length_method, between_method)
    
    # Helper function for rescaling the numeric data. 
    rescale_df_total = function(df, target) {
        s = mean(as.matrix(df), na.rm = TRUE) 
        scaled_num = df * (target / s)
        return(scaled_num)
    }
    
    # Choose common target.
    tot_ctrl_s = mean(as.matrix(IP_signal_sf_ctrl), na.rm = TRUE)
    tot_ctrl_n = mean(as.matrix(IP_noise_sf_ctrl), na.rm = TRUE)
    tot_trt_s = mean(as.matrix(IP_signal_sf_trt), na.rm = TRUE)
    tot_trt_n = mean(as.matrix(IP_noise_sf_trt), na.rm = TRUE)
    tot_in = mean(as.matrix(IN_sf), na.rm = TRUE)
    common_target = mean(c(tot_ctrl_n, tot_trt_n, tot_ctrl_s, tot_trt_s, tot_in))
    
    # Rescale.
    IP_signal_sf_ctrl_adj  = rescale_df_total(IP_signal_sf_ctrl, common_target)
    IP_noise_sf_ctrl_adj = rescale_df_total(IP_noise_sf_ctrl, common_target)
    IP_signal_sf_trt_adj = rescale_df_total(IP_signal_sf_trt, common_target)
    IP_noise_sf_trt_adj = rescale_df_total(IP_noise_sf_trt, common_target)
    IN_sf_adj = rescale_df_total(IN_sf, common_target)
    
    # Add the id back in (not needed for ctrl and trt since they will be combined).
    IN_sf_adj = rownames_to_column(IN_sf_adj, var = id_col)
    IP_signal_sf_ctrl_adj[[id_col]] = window_ctrl[[id_col]]
    IP_noise_sf_ctrl_adj[[id_col]] = background_ctrl[[id_col]]
    IP_signal_sf_trt_adj[[id_col]] = window_trt[[id_col]]
    IP_noise_sf_trt_adj[[id_col]] = background_trt[[id_col]]
    
    # Combine the ctrl and treatment stuff back into one group. 
    IP_sf_ctrl_adj = rbind(IP_signal_sf_ctrl_adj, IP_noise_sf_ctrl_adj)
    IP_sf_trt_adj = rbind(IP_signal_sf_trt_adj, IP_noise_sf_trt_adj)
    
    # Align the treatment to the control (ensures consistency). 
    ctrl_idx = rownames(IP_sf_ctrl_adj)
    IP_sf_trt_adj  = IP_sf_trt_adj[ctrl_idx,  , drop = FALSE]
    rownames(IP_sf_trt_adj)  = ctrl_idx
    
    # Remove ID from one of the treatment groups. 
    IP_sf_trt_adj[[id_col]] = NULL
    
    # Bind sample columns and attach merged gene_id.
    IP_sf_adj = cbind(IP_sf_ctrl_adj, IP_sf_trt_adj)

    # Add the window identifiers back in.
    IP_sf_adj = rownames_to_column(IP_sf_adj, var = "name")
    
    # Ensure that both id columns are characters so they can be joined on. 
    IP_sf_adj[[id_col]] = as.character(IP_sf_adj[[id_col]])
    IN_sf_adj[[id_col]] = as.character(IN_sf_adj[[id_col]])
    
    # Join all of the size factors together. 
    sf_adj = left_join(IP_sf_adj, IN_sf_adj, by = id_col)
    
    # Remove rows lacking complete normalization factors.
    sf_adj = sf_adj[complete.cases(sf_adj), , drop = FALSE]
    rownames(sf_adj) = NULL
    
    # Build the vector of column names for control IP replicates.
    ctrl_cols = setdiff(names(IP_sf_ctrl_adj), id_col)
    
    # Same for treatment IP replicates.
    trt_cols  = setdiff(names(IP_sf_trt_adj), id_col)

    # Extract the IP reps. 
    ctrl_ip_cols = get_reps_IP_only(gene_wise_IP, ctrl)
    trt_ip_cols  = get_reps_IP_only(gene_wise_IP, trt)
    
    # Create summed columns for calculating average background specific size factors. 
    gene_wise_IP[[paste0(ctrl, "_IP_sum")]] = rowSums(gene_wise_IP[ctrl_ip_cols], na.rm = TRUE)
    gene_wise_IP[[paste0(trt, "_IP_sum")]]  = rowSums(gene_wise_IP[trt_ip_cols],  na.rm = TRUE)
    
    # Size factors for average 
    ctrl_sum_col = paste0(ctrl, "_IP_sum")
    trt_sum_col  = paste0(trt, "_IP_sum")
    IP_backg_sf = get_sf_MOR(
          dplyr::select(gene_wise_IP, all_of(c("gene_id", ctrl_sum_col, trt_sum_col)))
        )
    
    # Combine the rep specific size factors with the average background size factors. 
    norm_factors = combine_ip_size_factors(sf_adj, IP_backg_sf, ctrl, trt)
    
    # Save the size factors. 
    write.table(norm_factors, output, row.names = FALSE, sep = "\t")
}
