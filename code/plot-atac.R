
plot_type_egc_scatter <- function(mcatac, cell_types, const_threshold = NULL, limits = c(-17, -12), pointsize = 2, plot_ablines = TRUE, plot_const_line = FALSE, const_line_value = -14.5, use_theme = TRUE, log_transform = TRUE, base_type = "Epiblast") {
    mc_egc <- as.matrix(mcatac@egc)
    if (log_transform){
        mc_egc <- log2(mc_egc + 1e-5)
    }

    # get cell type egc average for each peak
    ct_egc <- t(tgs_matrix_tapply(mc_egc, mcatac@metadata$cell_type, mean, na.rm = TRUE))

    ct_df <- ct_egc %>%
        as.data.frame() %>%
        rownames_to_column("cre")
    ct_df <- ct_df[, c("cre", cell_types)] %>%
        gather("cell_type", "egc", -cre) %>%
        left_join(ct_df %>% select(cre, !!sym(base_type)), by = "cre") %>%
        mutate(cell_type = factor(cell_type, levels = cell_types))

    if (!is.null(const_threshold)){
        const_v <- tibble(cre = rownames(ct_egc), const = matrixStats::rowMins(ct_egc) >= const_threshold) %>%
            select(cre, const) %>%
            deframe()
    } else {
        const_v <- mcatac@peaks$const
        if (is.null(const_v)){
            cli::cli_abort("{.field const_threshold} must be specified if {.field mcatac@peaks$const} is NULL")
        }
    }

    cts_colors <- mcatac@metadata %>%
        select(cell_type, cell_type_color) %>%
        distinct() %>%
        deframe()
    
    ct_df <- ct_df %>%
        mutate(x = !!sym(base_type)) 

    p <- ct_df %>%
        filter(!const_v[cre]) %>%        
        ggplot(aes(x = x, y = egc, color = cell_type)) +
        scattermore::geom_scattermore(data = ct_df %>% filter(const_v[cre]), pointsize = pointsize, color = "black") +
        scattermore::geom_scattermore(pointsize = pointsize) +
        scale_color_manual(values = cts_colors) +
        geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +        
        theme(legend.position = "none", aspect.ratio = 1) 
        
    if (!is.null(limits)) {
        p <- p +
            xlim(limits) +
            ylim(limits)
    }

    if (plot_ablines){
        p <- p +            
            geom_abline(slope = 1, intercept = -1, color = "black", linetype = "dashed") +
            geom_abline(slope = 1, intercept = 1, color = "black", linetype = "dashed") 
    }

    if (plot_const_line){
        p <- p +            
            geom_hline(yintercept = const_line_value, color = "red", linetype = "dashed")
    }

    if (length(cell_types) == 1) {
        p <- p + labs(x = base_type, y = cell_types)
    } else {
        p <- p +
            labs(x = base_type, y = "EGC") + 
            facet_wrap(~cell_type, ncol = 4, nrow = 4)
    }

    if (use_theme){
        p <- p + theme_classic() 
    }
        
    
    return(p)
}

plot_prob_cell_type_scatter <- function(ct_egc, cell_type1, cell_type2, similar_thresh = 0.1, low_thresh = 0.4){
    ct_df <- ct_egc %>%
        as.data.frame() %>%
        rownames_to_column("cre")
    
    plot_df <- ct_df[, c("cre", cell_type2)] %>%
        gather("cell_type", "egc", -cre) %>%
        left_join(ct_df %>% select(cre, !!sym(cell_type1)), by = "cre") %>%
        mutate(cell_type = factor(cell_type, levels = c(cell_type1, cell_type2))) %>%
        mutate(diff = egc - !!sym(cell_type1)) %>%
        mutate(type = case_when(
            abs(diff) <= similar_thresh ~ "similar",
            egc <= low_thresh & !!sym(cell_type1) <= low_thresh ~ "low",                        
            diff > similar_thresh ~ "higher",
            diff < -similar_thresh ~ "lower"
        ))

    p_prob_colored <- plot_df %>%
        ggplot(aes(x = !!sym(cell_type1), y = egc, color = type)) +
        scattermore::geom_scattermore(pointsize = 2) +
        ylab(cell_type2) +
        xlab(cell_type1) +
        scale_color_manual(values = c("low" = "orange", "similar" = "grey", "higher" = "red", "lower" = "#0077ff")) +
        guides(color = "none") +
        theme(aspect.ratio = 1) + 
        scale_y_continuous(breaks = c(0, 0.2, 0.5, 0.75, 1), labels = c("0", "0.2", "0.5", "0.75", "1"), limits = c(0, 1)) +
        scale_x_continuous(breaks = c(0, 0.2, 0.5, 0.75, 1), labels = c("0", "0.2", "0.5", "0.75", "1"), limits = c(0, 1))
        

    p_prob_colored
}

