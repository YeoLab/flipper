# Load necessary packages.
suppressPackageStartupMessages({
    library(glue)
    library(tidyverse)
    library(Matrix)
    library(parallel)
    library(tibble)
    library(dplyr)
    library(data.table)
})

# Load up the arguments passed from snakemake.
args = commandArgs(trailingOnly=TRUE)

# Divide arguments into proper names. 
sample_ctrl = args[1] 
sample_trt = args[2] 
skipper_dir = args[3]
anno_path = args[4]
substrate_control = args[5]
tool_dir = args[6]

# Load up functions from utils. 
source(glue("{tool_dir}/utils.R"))

# Load up count data. 
data = load_and_combine_skipper_count(skipper_dir, sample_ctrl, sample_trt, anno_path)
count_data = data$count_data
ctrl = data$ctrl
trt = data$trt
write.table(count_data, "intermediates/all_count_data.tsv", sep="\t", row.names = FALSE)

save(ctrl, trt, file = "intermediates/conditions.RData")

# Load up the window data. 
window_data = load_and_combine_window(skipper_dir, sample_ctrl, sample_trt)

# Get the individual count data from each rep for each window. 
window_counts_INw = make_INw(ctrl, trt, window_data, count_data)  

# Calculate genewise IN by summing over gene name in the count data. 
substrate_level_IN = get_substrate_level(count_data, ctrl, trt, substrate_control)$substrate_wise_IN

# Replace the window wise IN counts with the gene wise IN counts and save. 
replace_window_IN_count(window_counts_INw, substrate_level_IN, ctrl, trt, substrate_control)
