######################### Quality of life #########################

# Simple function that adjusts for Rs weird naming scheme (needed in case users start sample names with numbers). 
add_x_prefix = function(cond) {
    if (grepl("^[0-9]", cond)) {
        return(paste0("X", cond))
    } else {
        return(cond)
    }
}

# Simple helper function for finding the opposite of an intersection. 
outersect = function(x, y) {
    sort(union(setdiff(x, y), setdiff(y, x)))
}

# Various functions for extracting replicate names.
get_reps = function(df, condition) {
    grep(paste0("^", condition, "_(IN|IP)_"), names(df), value = TRUE)
}
get_reps_IN_only = function(df, condition) {
    grep(paste0("^", condition, "_IN_"), names(df), value = TRUE)
}
get_reps_IP_only = function(df, condition) {
    grep(paste0("^", condition, "_IP_"), names(df), value = TRUE)
}

######################### For preprocessing #########################

# automatically extracts out the condition names from the count data. 
extract_condition = function(count_data) {
    colnames_data = colnames(count_data)
    count_cols = colnames_data[8:length(colnames_data)]
    base_names = gsub("_IN_\\d+|_IP_\\d+", "", count_cols)
    cond_name = unique(base_names)
    return(cond_name)
}

load_and_combine_skipper_count = function(skipper_dir, sample_ctrl, sample_trt, anno_path) {
    
    # Load the complete count data. 
    count_ctrl = read.csv(paste0(skipper_dir, '/secondary_results/counts/genome/tables/',sample_ctrl, '.tsv.gz'),
                              sep = '\t')
    count_trt = read.csv(paste0(skipper_dir, '/secondary_results/counts/genome/tables/',sample_trt, '.tsv.gz'),
                              sep = '\t')

    # Extract out the ctrl and treatment names. 
    ctrl = extract_condition(count_ctrl)
    trt = extract_condition(count_trt)

    # Create a list of the ctrl columns of interest.  
    ctrl_cols  = get_reps(count_ctrl, ctrl)   
    select_col = c("chr","start","end","strand","gc","name", ctrl_cols)

    # Subset the count data to the columns of interest. 
    count_data = subset(count_ctrl, select=select_col)

    # Create a list of the treatment columns of interest. 
    trt_cols = get_reps(count_trt, trt)   
    
    # Add the count data 2 to count data 1
    count_data[trt_cols] = count_trt[trt_cols]

    # Add in length and gene id data. 
    annotations = read_tsv(anno_path)
    if (length(annotations$start) != length(count_data$start)) {
        stop(
            paste0(
                "Error: Mismatch in lengths.\n",
                "annotations: ", length(annotations$start), "\n",
                "count_data: ", length(count_data$start)
            )
        )
    }
    count_data$gene_id = annotations$gene_id
    count_data$base_feature_id = annotations$base_feature_id
    count_data$length = count_data$end - count_data$start
    
    return(list(count_data = count_data, ctrl = ctrl, trt = trt))
}

load_and_combine_window = function(skipper_dir, sample_ctrl, sample_trt) {

    # Load up the reproducible windows for each treatment group. 
    window_ctrl = read.csv(glue('{skipper_dir}/reproducible_enriched_windows/{sample_ctrl}.reproducible_enriched_windows.tsv.gz'),
                              sep = '\t')
    window_trt = read.csv(glue('{skipper_dir}/reproducible_enriched_windows/{sample_trt}.reproducible_enriched_windows.tsv.gz'),
                              sep = '\t')

    # Combine the ctrl and treatment windows. 
    window_data = rbind(window_ctrl,window_trt)

    # Subset to only the columns of interest. 
    window_data = subset(window_data, select = c(chr, start, end, name, score, strand, feature_type_top, feature_types, transcript_type_top, gene_name, gene_id))

    # Remove duplicated windows (no need to double count windows found in both). 
    window_data = window_data[!duplicated(window_data), ]
    return(window_data)
}

# Get the individual counts at every window. 
make_INw = function(ctrl, trt, window_data, count_data) {
    
    # Get the counts for all of the enriched window. 
    window_counts_INw = left_join(window_data, count_data)

    # Extract out control and treatment column names. 
    ctrl_cols = get_reps(window_counts_INw, ctrl)   
    trt_cols = get_reps(window_counts_INw, trt) 
    
    # establish a consistent column order for all datasets. 
    col_order = c('chr','start','end', 'name', 'score', 'strand','feature_type_top', 'feature_types',
                  'transcript_type_top','gene_name','gene_id','base_feature_id', ctrl_cols, trt_cols)
    window_counts_INw = window_counts_INw[, col_order]
    
    write.table(window_counts_INw, "intermediates/window_counts_INw.tsv", sep="\t", row.names = FALSE)
    return(window_counts_INw)
}

get_substrate_level = function(count_data, ctrl, trt, substrate_control) {
    
    # Establish what to sum the IN across.
    if (substrate_control == "INf") {
        id = "base_feature_id"
        other_id = "gene_id"
    } else if (substrate_control == "INg") {
        id = "gene_id"
        other_id = "base_feature_id"
    } else {
        stop("substrate_control must be either 'INf' or 'INg'.")
    }
    
    # Sum substrate-wise counts while keeping both IDs.
    substrate_wise_counts = count_data %>%
        group_by(.data[[id]]) %>%
        summarise(
            # Keep the grouping ID.
            !!id := dplyr::first(.data[[id]]),
            # Carry along the other ID.
            !!other_id := dplyr::first(.data[[other_id]]),
            # Summarize. 
            across(.cols = where(is.numeric) & !matches("^gc$"), .fns  = ~ sum(.x, na.rm = TRUE)),
            gc = mean(gc, na.rm = TRUE),
            .groups = "drop"
        )
    
    # Remove these broken columns. 
    substrate_wise_counts$start = NULL
    substrate_wise_counts$end = NULL

    # Subset to just the IN data. 
    trt_IN = get_reps_IN_only(substrate_wise_counts, trt)
    ctrl_IN =  get_reps_IN_only(substrate_wise_counts, ctrl)
    IN_subset = c("gene_id", "base_feature_id", "length", "gc", ctrl_IN, trt_IN)
    substrate_wise_IN = substrate_wise_counts[, IN_subset]

    # Subset to just the IP data.
    trt_IP = get_reps_IP_only(substrate_wise_counts, trt)
    ctrl_IP =  get_reps_IP_only(substrate_wise_counts, ctrl)
    IP_subset = c("gene_id", "base_feature_id", "length", "gc", ctrl_IP, trt_IP)
    substrate_wise_IP = substrate_wise_counts[, IP_subset]

    return(list(substrate_wise_IN = substrate_wise_IN, substrate_wise_IP = substrate_wise_IP))
}

# Replace the window wise IN counts (INw) with gene (INg) or feature (INf) wise IN counts. 
replace_window_IN_count = function(window_counts_INw, substrate_wise_IN, ctrl, trt, substrate_control) {

    # Establish what to sum the IN across.
    if (substrate_control == "INf") {
        id = "base_feature_id"
    } else if (substrate_control == "INg") {
        id = "gene_id"
    } else {
        stop("substrate_control must be either 'INf' or 'INg'.")
    }

    # Define the columns to keep for the final dataset. 
    trt_cols = get_reps(window_counts_INw, trt)
    ctrl_cols = get_reps(window_counts_INw, ctrl)
    col_order = c('chr','start','end', 'name', 'score', 'strand', 'feature_type_top', 'feature_types', 
        'transcript_type_top','gene_name', 'gene_id', 'base_feature_id', trt_cols, ctrl_cols)

    # Get a list of the INw columns to drop. 
    trt_cols = get_reps_IN_only(window_counts_INw, trt)
    ctrl_cols = get_reps_IN_only(window_counts_INw, ctrl)
    drops = c(trt_cols, ctrl_cols)
    
    # Replace the window-wise IN counts with the gene/feature-wise IN. 
    temp = window_counts_INw[ , !(names(window_counts_INw) %in% drops)]

    # Subset to just the ID columns. 
    substrate_wise_IN_clean = substrate_wise_IN %>%
    select(-any_of(setdiff(intersect(names(temp), names(substrate_wise_IN)), id)))

    # Left join to add the new IN data back into the window counts. 
    window_counts_INsub = left_join(temp, substrate_wise_IN_clean, by = id)

    # Fix column order. 
    window_counts_INsub = window_counts_INsub[, col_order]
    
    write.table(window_counts_INsub, "intermediates/window_counts_INsub.tsv", sep="\t", row.names = FALSE)
    return(list(window_counts_INsub = window_counts_INsub))
}

######################### For size factors #########################

# Loads up IP windows and splits them into background windows and significant binding windows. 
get_window_counts = function(count_data, skipper_dir, sample_cond, cond) {
    windows = read.csv(paste0(skipper_dir, '/reproducible_enriched_windows/',sample_cond,'.reproducible_enriched_windows.tsv.gz'),
                              sep = '\t')
    trt_cols = get_reps_IP_only(count_data, cond)   
    IP_cols = c(trt_cols, "gene_id", "base_feature_id", "length", "gc", "name")
    window_counts = dplyr::select(left_join(windows, count_data), all_of(IP_cols))
    background_counts = dplyr::select(anti_join(count_data, windows), all_of(IP_cols))
    return(list(window_counts = window_counts, background_counts = background_counts))
}

# Rescale a size-factor vector so its mean equals 1. 
rescale_to_mean = function(x, target = 1) {
    x * (target / mean(x, na.rm = TRUE))
}

# Rescale a numeric data frame so its total sum equals a target sum. 
rescale_df_total = function(df, target) {
    s = mean(as.matrix(df), na.rm = TRUE) 
    scaled_num = df * (target / s)
    return(scaled_num)
}

combine_ip_size_factors = function(norm_factors, IP_backg_sf, ctrl, trt, ip_pattern = "_IP_") {
    
    # Helper for geometric mean.
    GM = function(x) exp(mean(log(x)))
    
    # Get the sample/replicate names.
    nms = names(norm_factors)
    
    # Find IP replicate columns for ctrl / trt.
    ctrl_ip_cols = grep(paste0("^", ctrl, ip_pattern), nms, value = TRUE)
    trt_ip_cols  = grep(paste0("^", trt,  ip_pattern), nms, value = TRUE)
    
    # Pull summed IP size factors.
    ctrl_sum_name = paste0(ctrl, ip_pattern, "sum")
    trt_sum_name  = paste0(trt,  ip_pattern, "sum")
    
    ctrl_sum_sf = IP_backg_sf[ctrl_sum_name]
    trt_sum_sf  = IP_backg_sf[trt_sum_name]
    
    # Target ratio from summed libraries.
    R = as.numeric(ctrl_sum_sf / trt_sum_sf)

    if (is.null(dim(norm_factors))) {
        d = GM(norm_factors[ctrl_ip_cols])
        o = GM(norm_factors[trt_ip_cols])

        a = sqrt(R * (o / d))
        b = 1 / a

        scaled = norm_factors
        scaled[ctrl_ip_cols] = a * scaled[ctrl_ip_cols]
        scaled[trt_ip_cols]  = b * scaled[trt_ip_cols]        

    } else {
        d = as.numeric(GM(colMeans(norm_factors[, ctrl_ip_cols, drop = FALSE])))
        o = as.numeric(GM(colMeans(norm_factors[, trt_ip_cols,  drop = FALSE])))

        a = sqrt(R * (o / d))
        b = 1 / a

        scaled = norm_factors
        scaled[, ctrl_ip_cols] = a * scaled[, ctrl_ip_cols, drop = FALSE]
        scaled[, trt_ip_cols]  = b * scaled[, trt_ip_cols,  drop = FALSE]        
        }
    return(scaled)
}

get_sf_MOR = function(data) {
    
    # Pick only count columns. 
    count_cols = setdiff(names(data),
                          c("gene_name","length","gc","name","gene_id"))
    
    # Coerce to numeric matrix.
    counts_df = as.data.frame(data[, count_cols, drop = FALSE])
    counts_df[] = lapply(counts_df, function(x) as.numeric(as.character(x)))
    counts = as.matrix(counts_df)
    
    # Set gene rownames. 
    rownames(counts) = as.character(data$gene_id)
    
    # Remove rows with any NA. 
    ok = complete.cases(counts)
    counts = counts[ok, , drop = FALSE]
    storage.mode(counts) = "numeric"  # ensure numeric
    
    # Get size factors. 
    sf = DESeq2::estimateSizeFactorsForMatrix(counts, type = "poscounts")
    return(sf)
}

# Drop all rows with zero counts (needed for EDAseq). 
drop_0_rows = function(data) {
    just_counts = data[, setdiff(names(data), c("gene_id", "base_feature_id", "name", "gene_name", "length", "gc"))]
    keep = rowSums(as.matrix(just_counts)) > 0
    data = data[keep, , drop = FALSE]
    rownames(data) = NULL
    return(data)
}

# If there are 0 counts in a row from any of these datasets, drop them in all datasets. 
drop_0_rows_joint_universe = function(window_ctrl, background_ctrl, window_trt, background_trt) {
    
    # Grab the metadata columns (to be ignored when counting)
    metadata_cols = c("gene_id", "base_feature_id", "name", "gene_name", "length", "gc")

    # Helper function for finding window name of windows with 0 counts. 
    get_nonzero_names = function(data) {
        count_cols = setdiff(names(data), metadata_cols)
        counts = data[, count_cols, drop = FALSE]
        return(as.character(data$name[rowSums(as.matrix(counts), na.rm = TRUE) > 0]))
    }

    # Find all of the names that have at least 1 count in each dataset. 
    keep_names = Reduce(union,list(get_nonzero_names(window_ctrl),get_nonzero_names(background_ctrl),
                                    get_nonzero_names(window_trt),get_nonzero_names(background_trt)))

    return(list(window_ctrl = window_ctrl[as.character(window_ctrl$name) %in% keep_names, , drop = FALSE],
                background_ctrl = background_ctrl[as.character(background_ctrl$name) %in% keep_names, , drop = FALSE],
                window_trt = window_trt[as.character(window_trt$name) %in% keep_names, , drop = FALSE],
                background_trt = background_trt[as.character(background_trt$name) %in% keep_names, , drop = FALSE]))
}

get_sf_EDA = function(conds, data, identifier, reps, gc_method, length_method, between_method) {
   
    # Find the count columns. 
    count_cols = unlist(lapply(conds, function(cond) get_reps(data, cond)))

    # Subset to just the identifier and the counts. 
    identLevelData = data[, c(identifier, count_cols)]

    # Subset to just the features of interest and the counts. 
    feature = data[, c(identifier, "length", "gc")]

    # Reset to default rownames (can cause problems for column_to_rownames).
    rownames(identLevelData) = NULL
    rownames(feature) = NULL

    # Add identifier as the rownames. 
    identLevelData = identLevelData %>% column_to_rownames(var = identifier)
    feature = feature %>% column_to_rownames(var = identifier)

    # Run EDAseq. 
    EDAseq_results = newSeqExpressionSet(
        counts = as.matrix(identLevelData),
        featureData = feature,
        phenoData = data.frame(
            conditions = factor(rep(conds, each = reps)),
            row.names = colnames(identLevelData)
        )
    )

    # Run each of the different normalization types for EDAseq.
    dataOffset = withinLaneNormalization(EDAseq_results, "gc", which=gc_method, offset=TRUE)
    if (identifier %in% c("gene_id", "base_feature_id")) {
        dataOffset = withinLaneNormalization(dataOffset, "length", which=length_method, offset=TRUE)
    } # We do not want to run length normalization on windows. 
    dataOffset = betweenLaneNormalization(dataOffset, which=between_method, offset=TRUE)

    # Perform the adjustments reccomended by EDAseq. 
    offexp = exp(-1 * as.matrix(offst(dataOffset)))
    off_adj = offexp / exp(rowMeans(log(offexp)))

    return(data.frame(off_adj))
}

######################### For differential analysis #########################
                          
remove_na_rows = function(sf, counts) {

    # Find complete rows in each dataframe.
    complete_sf   = complete.cases(sf)
    complete_counts = complete.cases(counts)
    
    # Keep rows that are complete in both.
    keep = complete_sf & complete_counts
    
    return(list(sf = sf[keep, , drop = FALSE], counts = counts[keep, , drop = FALSE]))
}

process_differential = function(window_counts, sf, ctrl, trt, select_cols, pval_threshold, log2fc_threshold){

    if (is.null(dim(sf))) {
        wc_sub = tidyr::drop_na(window_counts)
        wc_mat = as.matrix(dplyr::select(wc_sub, all_of(select_cols)))

        # Ensure proper ordering. 
        sf = sf[select_cols]
        wc_mat = wc_mat[, select_cols, drop = FALSE]
    } else {
        # Decide which windows to keep (shared between the two tables).
        keep_windows = intersect(sf$name, window_counts$name)
        
        # Reorder both tables to the same window order, then subset to the same sample columns.
        sf_sub = sf[match(keep_windows, sf$name), colnames(sf), drop = FALSE]
        wc_sub = window_counts[match(keep_windows, window_counts$name), colnames(window_counts), drop = FALSE]
        
        # Convert to matrices.
        sf_mat = as.matrix(dplyr::select(sf_sub, all_of(select_cols)))
        wc_mat = as.matrix(dplyr::select(wc_sub, all_of(select_cols)))
        
        # Remove rownames
        rownames(sf_mat) = NULL
        rownames(wc_mat) = NULL

        # Ensure proper ordering. 
        sf_mat = sf_mat[, select_cols, drop = FALSE]
        wc_mat = wc_mat[, select_cols, drop = FALSE]
    }

    # Define the model and design for deseq. 
    sample_names = colnames(wc_mat)
    model = data.frame(
        Reps = paste0("Rep", sub(".*_(IN|IP)_", "", sample_names)),
        IP   = ifelse(grepl("_IN_", sample_names), "Input",
                    ifelse(grepl("_IP_", sample_names), "IP", NA)),
        Trt  = ifelse(grepl(paste0("^", trt, "_"),  sample_names), "Trt",
                    ifelse(grepl(paste0("^", ctrl, "_"), sample_names), "Ctrl", NA)),
        row.names = sample_names,
        stringsAsFactors = FALSE
    )
    design = ~Reps + IP*Trt

    # Create the deseq object. 
    dds = DESeqDataSetFromMatrix(countData = wc_mat,
                                colData = model,
                                design = design)

    # If sf is a vector, treat it as size factors. If a matrix, as normalization factors.
    if (is.null(dim(sf))) {
        sizeFactors(dds) = sf
    } else {
        normalizationFactors(dds) = sf_mat # Note that the offset alterations reccomended by EDAseq were already performed in calc_norm. 
    }

    # Run deseq.
    dds = DESeq2::estimateDispersions(dds, fitType = "local")
    dds = DESeq2::nbinomWaldTest(dds)

    # Extract deseq results into a cleaned table.
    res = DESeq2::results(dds, name = "IPIP.TrtTrt") 
    dat = data.frame(padj = res$padj, pval = res$pvalue, stat = res$stat, log2FoldChange = res$log2FoldChange)

    # Extract deseq results into a cleaned table. 
    result_check = cbind(wc_sub, dat)
    
    # Sort be the adjusted p-value. 
    result_check = result_check %>% arrange(padj)

    # Extract out the significant results. 
    sigs = dplyr::filter(result_check, padj < pval_threshold, abs(log2FoldChange) > log2fc_threshold)

    # Return the final results. 
    return(list(dds = dds, window_counts = window_counts, wc_mat = wc_mat, result_check = result_check, sigs = sigs))
}

# Function to compute Fisher’s combined p-value.
combine_p_fisher = function(pvals) {
    stat = -2 * sum(log(pvals))
    pchisq(stat, df = 2 * length(pvals), lower.tail = FALSE)
}

######################### For plotting #########################

volcano_plot = function(data, padj_threshold = 0.05, log2fc_threshold = 1.0, top_n = 10,
                        colors = NULL, plot_width = 1, inside_legend = FALSE) {

    # Make sure data is the proper format. 
    data = as.data.frame(data)

    # Create empty plot in the event that there are no significant hits. 
    if (nrow(data) < 1) {
        p = ggplot() + theme_void() +
            labs(title = "Error", subtitle = "No data available (Treatment/Control)") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
        return(p)
    }
    
    # Compute transformed p-values and add new significance column. 
    data$neg_log10_padj = -log10(pmax(data$padj, 0))
    data$significance = "Not Significant"
    data$significance[data$padj < padj_threshold & data$log2FoldChange > log2fc_threshold] = "Up"
    data$significance[data$padj < padj_threshold & data$log2FoldChange < -1*log2fc_threshold] = "Down"
    data$significance[data$padj < padj_threshold & abs(data$log2FoldChange) <= log2fc_threshold] = "Pass p-value cutoff"
    
    # Top genes for labeling. 
    top_genes = data[data$significance != "Not Significant", ]
    top_genes = top_genes[order(top_genes$padj), ]
    top_genes = head(top_genes, top_n)
    top_genes$label = top_genes$gene_name
    
    # Symmetric x-limit.
    max_l2fc = max(abs(data$log2FoldChange), na.rm = TRUE)

    # Create Y-limit designed to contain a majority of the data without 
    # Being overly affected by outliers. 
    vals = data$neg_log10_padj[is.finite(data$neg_log10_padj)]
    if (length(vals) > 0) {
        q1 = quantile(vals, 0.05, na.rm = TRUE)
        q3 = quantile(vals, 0.95, na.rm = TRUE)
        iqr = q3 - q1
        fence = q3 + 8 * iqr
        capped = vals[vals <= fence]
        if (length(capped) == 0) {
        ymax = max(vals)
        } else {
        ymax = max(capped) * 1.1
        }
        ymax = max(ymax, 3)
    } else {
      ymax = 3
    }
    
    # Build plot.
    p = ggplot(data, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
        geom_point(alpha = 0.5, size = 2) +
        scale_color_manual(values = c("Not Significant" = "gray70", "Up" = "red",
                                    "Down" = "blue", "Pass p-value cutoff" = "black"),
                         name = "Significance (Treatment/Control)" ) +
        geom_vline(xintercept = c(-1*log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", linewidth = 0.5) +
        labs(x = expression("Log"[2]*" fold change"), y = expression("-Log"[10]*"(adjusted p-value)")) +
        coord_cartesian(xlim = c(-max_l2fc, max_l2fc), ylim = c(0, ymax)) +
        theme_minimal() +
        theme(panel.grid.major = element_line(color = "gray90", linewidth = 0.3), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA),
            plot.title = element_text(hjust = 0.5, face = "bold"), axis.title = element_text(size = 12),
            axis.text = element_text(size = 10), legend.position = "right", legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 8), legend.key.size = unit(0.5, "cm"))
  
    return(p)
}

feature_type_bar_plot = function(data, padj_threshold = 0.05, log2fc_threshold = 1.0) {
    
    # Make sure data is the proper format. 
    df = as.data.frame(data)

    # Add a new significance column. 
    df$significance = NA
    df$significance[df$padj < padj_threshold & df$stat >  log2fc_threshold]  = "Up"
    df$significance[df$padj < padj_threshold & df$stat < -1*log2fc_threshold] = "Down"
    sig = df %>% filter(significance %in% c("Up", "Down"))
    
    # Create empty plot in the event that there are no significant hits. 
    if (nrow(sig) < 1) {
        return(
          ggplot() +
            theme_void() +
            labs(title = "Distribution of Feature Types", subtitle = "No data available") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
        )
    }
    
    # Calculate counts per feature.
    counts = sig %>%
        dplyr::count(feature_type_top, significance) %>%
        tidyr::pivot_wider(names_from = significance, values_from = n, values_fill = 0)
    
    # Ensure both Up and Down columns exist. 
    if (!"Up" %in% names(counts)) {
      counts$Up = 0
    }
    if (!"Down" %in% names(counts)) {
      counts$Down = 0
    }
    
    counts = counts %>%
      mutate(Down = -Down,
             total = abs(Down) + Up) %>%
      arrange(desc(total))
    
    # Prepare for plotting.
    plot_data = counts %>%
      select(feature_type_top, Up, Down) %>%
      pivot_longer(cols = c(Up, Down), names_to = "Regulation", values_to = "Count") %>%
      mutate(feature_type_top = factor(feature_type_top, levels = counts$feature_type_top))
    
    # Compute a small offset for label nudging.
    max_c = max(abs(plot_data$Count))
    offset = max_c * 0.02
    
    # Build plot.
    p = ggplot(plot_data, aes(x = feature_type_top, y = Count, fill = Regulation)) +
        geom_col(width = 0.8) +
        coord_flip(clip = "off") +
        scale_fill_manual(values = c("Down" = "blue", "Up" = "red")) +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        # Only label non-zero values.
        geom_text(data = subset(plot_data, Count > 0), aes(y = Count + offset, label = Count), hjust = 0, size = 3) +
        geom_text(data = subset(plot_data, Count < 0), aes(y = Count - offset, label = abs(Count)), hjust = 1, size = 3) +
        labs(title = "Distribution of Feature Types (Treatment/Control)", x = "Feature Type",
             y = "Number of Binding Windows", fill = "Regulation") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.y = element_text(size = 8, face = "bold"),
              axis.title = element_text(size = 12), panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(color = "gray90", size = 0.3), legend.position = "top",
              plot.margin = margin(t = 10, r = 40, b = 10, l = 10))
    return(p)
}

plot_top_genes = function(df, score_col, lfc_col, n = 20, score_label = NULL) {
  
    # Add unique plotting labels where gene symbols are duplicated.
    df = df %>% mutate(gene_label = if_else(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE), 
                                          paste0(gene_name, " (", gene_id, ")"), gene_name))
    
    # Find the top upregulated and downregulated genes. 
    down = df %>% filter(!!sym(lfc_col) < 0) %>% arrange(!!sym(score_col)) %>% slice_head(n = n)
    up = df %>% filter(!!sym(lfc_col) > 0) %>% arrange(!!sym(score_col)) %>% slice_head(n = n)

    # Combine top up and down into "top"
    top_genes = bind_rows(down, up) %>% arrange(desc(!!sym(lfc_col))) %>%
    mutate(gene_label = factor(gene_label, levels = rev(unique(gene_label))), score_log10 = -log10(!!sym(score_col)),
           score_log10 = ifelse(is.infinite(score_log10), 300, score_log10))

    # Find the maximum absolute value for y-limits. 
    max_abs = max(abs(top_genes[[lfc_col]]), na.rm = TRUE)

    # Find the minimum (closest to 0) of the up and downregulated so that we can draw a line between them. 
    neg_row = top_genes %>% filter(!!sym(lfc_col) < 0) %>% slice_max(!!sym(lfc_col), n = 1)
    pos_row = top_genes %>% filter(!!sym(lfc_col) > 0) %>% slice_min(!!sym(lfc_col), n = 1)

    # Used to prevent outliers from messing with the colour scale limit. 
    raw_max = max(top_genes$score_log10, na.rm = TRUE)
    p9 = quantile(top_genes$score_log10, 0.9, na.rm = TRUE)
    col_max = min(raw_max, p9 * 1.5)
    col_max = max(col_max, 1)

    # Create the plot. 
    p = ggplot(top_genes, aes(x = !!sym(lfc_col), y = gene_label)) + 
        geom_point(aes(size = n_sites, color = score_log10)) +
        scale_color_gradient(low = "blue", high = "red", limits = c(0, col_max), oob = scales::squish, 
                             name = if (is.null(score_label)) paste0("-log10(", score_col, ")") else score_label) +
        scale_size_continuous(name = "Number of sites", range = c(3, 8)) +
        scale_x_continuous(limits = c(-max_abs, max_abs)) +
        geom_vline(xintercept = 0, color = "red") +
        {if (nrow(neg_row) > 0 && nrow(pos_row) > 0)
          geom_segment(
            data = NULL,
            aes(
              x = neg_row[[lfc_col]], y = neg_row$gene_label,
              xend = pos_row[[lfc_col]], yend = pos_row$gene_label
            ),
            inherit.aes = FALSE,
            linetype = "dotted", color = "grey50"
          )} +
        theme_minimal(base_size = 17) +
        theme(axis.title.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
        labs(x = ifelse(lfc_col == "mean_logFC", "Mean of log2 Fold Change", "log2 Fold Change at Highest Confidence Site"))
  
    return(p)
}
