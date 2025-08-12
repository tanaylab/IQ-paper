# Helper functions for Figure 5 analysis

# Function to compare any two models with paired statistical tests
compare_models_multiple <- function(data, models,
                                    title = NULL,
                                    model_names_map = NULL,
                                    organism_colors = c("Human" = "red", "Mouse" = "#030366"),
                                    class_colors = c("IQ Ensemble" = "#3e6622", "Ensemble" = "#7575d4", "CNN" = "#CE3D32FF", "Transformer" = "#F0E685FF", "IQ" = "#749B58FF"),
                                    ylim_range = NULL,
                                    show_cell_types = TRUE,
                                    hide_ns_bars = TRUE,
                                    p_adjust_method = "BH",
                                    test_method = "wilcox",
                                    show_exact_pval = TRUE,
                                     subtitle = NULL, 
                                    show_connected_lines = TRUE,
                                    linetype = "solid",
                                    line_color = "black",
                                    line_alpha = 0.6,
                                    line_width = 0.8,
                                    point_size = 2.5,
                                    show_signif = TRUE,
                                    test_all_pairs = FALSE,
                                    explicit_comparisons = NULL, 
                                     return_pairwise_results = FALSE,
                                     conf_level = 0.95) { 
    # Validate input models
    if (!all(models %in% unique(data$model))) {
        missing_models <- models[!models %in% unique(data$model)]
        stop(paste(
            "Some models not found in data:",
            paste(missing_models, collapse = ", "),
            "\nAvailable models:",
            paste(unique(data$model), collapse = ", ")
        ))
    }

    if (length(models) < 2) {
        stop("At least two models must be provided for comparison")
    }

    # Validate test method
    valid_methods <- c("wilcox", "t.test", "permutation", "bootstrap")
    if (!(test_method %in% valid_methods)) {
        stop(paste("Invalid test method. Choose from:", paste(valid_methods, collapse = ", ")))
    }

    # Filter data to include only the specified models
    model_data <- data %>%
        filter(model %in% models)

    # Ensure model factor levels are in the order provided
    model_data$model <- factor(model_data$model, levels = models)

    # Apply model name mapping if provided
    if (!is.null(model_names_map)) {
        # Create a mapping just for the models being compared
        display_names <- model_names_map[models]
        model_data$model_display <- factor(model_names_map[as.character(model_data$model)],
            levels = display_names
        )
    } else {
        model_data$model_display <- model_data$model
    }

    # Reshape data for statistical tests
    wide_data <- model_data %>%
        select(cell_type, organism, model, r2) %>%
        pivot_wider(
            id_cols = c(cell_type, organism),
            names_from = model,
            values_from = r2
        )

    # Set default title if not provided
    if (is.null(title)) {
        title <- paste("Comparison of", paste(model_names_map[models], collapse = ", "), "Models")
    }

    # Calculate max y for annotation placement
    max_y <- max(model_data$r2, na.rm = TRUE)

    # Set default y-limit if not provided
    if (is.null(ylim_range)) {
        ylim_range <- c(0, max_y * 1.15) # Slightly higher to accommodate more annotations
    }

    # Function to perform permutation test
    perform_permutation_test <- function(x, y, n_permutations = 1000) {
        observed_diff <- mean(y - x, na.rm = TRUE)

        # Store permutation differences
        perm_diffs <- numeric(n_permutations)

        # Perform permutations
        for (i in 1:n_permutations) {
            # Randomly flip signs of differences
            signs <- sample(c(-1, 1), length(x), replace = TRUE)
            perm_diffs[i] <- mean(signs * (y - x), na.rm = TRUE)
        }

        # Calculate p-value (two-sided)
        p_value <- mean(abs(perm_diffs) >= abs(observed_diff))

        return(list(
            p_value = p_value,
            observed_diff = observed_diff
        ))
    }

    # Function to perform bootstrap test
    perform_bootstrap_test <- function(x, y, n_bootstraps = 1000) {
        observed_diff <- mean(y - x, na.rm = TRUE)

        # Calculate bootstrap distribution of differences
        boot_diffs <- numeric(n_bootstraps)
        n <- length(x)

        for (i in 1:n_bootstraps) {
            # Bootstrap sampling with replacement
            indices <- sample(1:n, n, replace = TRUE)
            boot_x <- x[indices]
            boot_y <- y[indices]
            boot_diffs[i] <- mean(boot_y - boot_x, na.rm = TRUE)
        }

        # Calculate 95% confidence interval
        ci <- quantile(boot_diffs, c(0.025, 0.975))

        # Calculate p-value based on whether CI includes 0
        p_value <- ifelse(ci[1] <= 0 && ci[2] >= 0, 1, 0.001)

        return(list(
            p_value = p_value,
            observed_diff = observed_diff,
            ci_lower = ci[1],
            ci_upper = ci[2]
        ))
    }

    # Perform pairwise statistical tests
    pairwise_results <- list()
    pair_counter <- 1
    
    # Determine which pairs to test
    if (!is.null(explicit_comparisons)) {
        # Use explicitly specified comparisons
        # Expected format: list(c("model1", "model2"), c("model1", "model3"), ...)
        if (!is.list(explicit_comparisons) || 
            !all(sapply(explicit_comparisons, function(pair) {
                is.vector(pair) && length(pair) == 2 && all(pair %in% models)
            }))) {
            stop("explicit_comparisons must be a list of 2-element vectors, each containing valid model names")
        }
        
        model_pairs <- do.call(rbind, lapply(explicit_comparisons, function(pair) {
            data.frame(model_a = pair[1], model_b = pair[2], stringsAsFactors = FALSE)
        }))
    } else if (test_all_pairs) {
        # Test all possible pairs
        model_pairs <- expand.grid(
            model_a = models,
            model_b = models,
            stringsAsFactors = FALSE
        ) %>%
        filter(model_a != model_b) %>%
        # Only test each pair once (A vs B, not also B vs A)
        filter(match(model_a, models) < match(model_b, models))
    } else {
        # Test only adjacent pairs in the order provided
        model_pairs <- data.frame(
            model_a = models[1:(length(models) - 1)],
            model_b = models[2:length(models)],
            stringsAsFactors = FALSE
        )
    }
    
    # Perform tests for selected pairs
    for (i in 1:nrow(model_pairs)) {
        model_a <- model_pairs$model_a[i]
        model_b <- model_pairs$model_b[i]
        
        # Get positions in the original models vector
        pos_a <- match(model_a, models)
        pos_b <- match(model_b, models)

        x <- wide_data[[model_a]]
        y <- wide_data[[model_b]]

        # Skip if all values are NA
        if (all(is.na(x)) || all(is.na(y))) {
            next
        }

        # Remove paired observations where either is NA
        valid_indices <- !is.na(x) & !is.na(y)
        x_valid <- x[valid_indices]
        y_valid <- y[valid_indices]

        # Skip if no valid pairs remain
        if (length(x_valid) == 0) {
            next
        }

        # Perform selected statistical test
        if (test_method == "wilcox") {
            # Wilcoxon signed-rank test (non-parametric, suitable for R² values)
            test_result <- wilcox.test(x_valid, y_valid, paired = TRUE, exact = FALSE, conf.int = TRUE, conf.level = conf_level)
            p_value <- test_result$p.value
            test_name <- "Wilcoxon signed-rank test"
            statistic_name <- names(test_result$statistic)
            statistic_value <- as.numeric(test_result$statistic)
            # Approximate z from two-sided p-value; sign by median difference
            median_diff <- stats::median(y_valid - x_valid, na.rm = TRUE)
            z_value <- stats::qnorm(p_value / 2, lower.tail = FALSE)
            z_value <- z_value * sign(median_diff)
            effect_size_r <- z_value / sqrt(length(x_valid))
            estimate <- if (!is.null(test_result$estimate)) as.numeric(test_result$estimate) else NA_real_
            conf_int_low <- if (!is.null(test_result$conf.int)) as.numeric(test_result$conf.int[1]) else NA_real_
            conf_int_high <- if (!is.null(test_result$conf.int)) as.numeric(test_result$conf.int[2]) else NA_real_
            degrees_freedom <- NA_real_
        } else if (test_method == "t.test") {
            # Paired t-test (parametric)
            test_result <- t.test(x_valid, y_valid, paired = TRUE, conf.level = conf_level)
            p_value <- test_result$p.value
            test_name <- "Paired t-test"
            statistic_name <- names(test_result$statistic)
            statistic_value <- as.numeric(test_result$statistic)
            degrees_freedom <- as.numeric(test_result$parameter)
            # Cohen's dz for paired samples
            diffs <- y_valid - x_valid
            effect_size_r <- NA_real_
            # Convert t to r as an effect size alternative: r = t / sqrt(t^2 + df)
            r_from_t <- statistic_value / sqrt(statistic_value^2 + degrees_freedom)
            # Keep both; store r_from_t in effect_size_r
            effect_size_r <- r_from_t
            estimate <- if (!is.null(test_result$estimate)) as.numeric(test_result$estimate) else mean(diffs, na.rm = TRUE)
            conf_int_low <- if (!is.null(test_result$conf.int)) as.numeric(test_result$conf.int[1]) else NA_real_
            conf_int_high <- if (!is.null(test_result$conf.int)) as.numeric(test_result$conf.int[2]) else NA_real_
            z_value <- NA_real_
        } else if (test_method == "permutation") {
            # Permutation test (robust, distribution-free)
            test_result <- perform_permutation_test(x_valid, y_valid)
            p_value <- test_result$p_value
            test_name <- "Permutation test"
            statistic_name <- "mean_diff"
            statistic_value <- as.numeric(test_result$observed_diff)
            degrees_freedom <- NA_real_
            estimate <- statistic_value
            conf_int_low <- NA_real_
            conf_int_high <- NA_real_
            z_value <- NA_real_
            effect_size_r <- NA_real_
        } else if (test_method == "bootstrap") {
            # Bootstrap test (robust confidence intervals)
            test_result <- perform_bootstrap_test(x_valid, y_valid)
            p_value <- test_result$p_value
            test_name <- "Bootstrap test"
            statistic_name <- "mean_diff"
            statistic_value <- as.numeric(test_result$observed_diff)
            degrees_freedom <- NA_real_
            estimate <- statistic_value
            conf_int_low <- as.numeric(test_result$ci_lower)
            conf_int_high <- as.numeric(test_result$ci_upper)
            z_value <- NA_real_
            effect_size_r <- NA_real_
        }

        # Store results
        pairwise_results[[paste(model_a, model_b, sep = "_vs_")]] <- list(
            model_a = model_a,
            model_b = model_b,
            p_value = p_value,
            position_a = pos_a,
            position_b = pos_b,
            n_samples = length(x_valid),
            test_name = test_name,
            statistic_name = statistic_name,
            statistic_value = statistic_value,
            degrees_freedom = degrees_freedom,
            z_value = z_value,
            effect_size_r = effect_size_r,
            estimate = estimate,
            conf_int_low = conf_int_low,
            conf_int_high = conf_int_high,
            conf_level = conf_level
        )
        
        pair_counter <- pair_counter + 1
    }

    # Extract p-values for adjustment
    p_values <- sapply(pairwise_results, function(x) x$p_value)

    # Adjust p-values for multiple comparisons
    adjusted_p_values <- p.adjust(p_values, method = p_adjust_method)

    # Update pairwise results with adjusted p-values
    for (i in 1:length(pairwise_results)) {
        pairwise_results[[i]]$adjusted_p_value <- adjusted_p_values[i]

        # Determine significance
        p_val <- pairwise_results[[i]]$adjusted_p_value
        pairwise_results[[i]]$significance <- if (p_val < 0.001) {
            "***"
        } else if (p_val < 0.01) {
            "**"
        } else if (p_val < 0.05) {
            "*"
        } else {
            "ns"
        }

        pairwise_results[[i]]$is_significant <- p_val < 0.05

        # Calculate annotation y-position with some spacing between annotations
        spacing_factor <- 0.05 * (i - 1)
        pairwise_results[[i]]$annotation_y <- max_y * (1.05 + spacing_factor)
    }

    if (return_pairwise_results) {
        # Build a rich results data frame including statistic, df (or n), p-values, effect size, and CIs
        p_res <- purrr::imap_dfr(
            pairwise_results,
            ~ {
                model_key <- .y
                res <- .x
                # Degrees of freedom for reporting: for t-test it's df; for Wilcoxon we use n-1 by convention, though Wilcoxon has no df.
                df_report <- if (!is.na(res$degrees_freedom)) {
                    res$degrees_freedom
                } else if (!is.na(res$n_samples)) {
                    max(res$n_samples - 1, 1)
                } else {
                    NA_real_
                }
                effect_label <- if (grepl("Wilcoxon|t-test", res$test_name)) "r" else "r"
                ci_percent <- round(res$conf_level * 100)
                # Format numbers
                fmt <- function(x) ifelse(is.na(x), NA_character_, formatC(x, digits = 3, format = "f"))
                stat_label <- ifelse(is.null(res$statistic_name), "stat", res$statistic_name)
                # Build APA-like string with df value (for non-parametric tests, df is reported as n-1 by convention)
                apa <- paste0(
                    stat_label, "(", ifelse(is.na(df_report), "NA", fmt(df_report)), ") = ", fmt(res$statistic_value),
                    ", p = ", fmt(res$p_value),
                    ", ", effect_label, " = ", fmt(res$effect_size_r),
                    ", ", ci_percent, "% Confidence Intervals = [", fmt(res$conf_int_low), ", ", fmt(res$conf_int_high), "]"
                )
                # Strict APA-like string using adjusted p-value as the primary p if multiple comparisons are present
                apa_strict <- paste0(
                    stat_label, "(", ifelse(is.na(df_report), "NA", fmt(df_report)), ") = ", fmt(res$statistic_value),
                    ", p = ", fmt(res$adjusted_p_value),
                    ", ", effect_label, " = ", fmt(res$effect_size_r),
                    ", ", ci_percent, "% Confidence Intervals = [", fmt(res$conf_int_low), ", ", fmt(res$conf_int_high), "]"
                )
                tibble(
                    model = model_key,
                    model_a = sub("_vs_.*$", "", model_key),
                    model_b = sub("^.*_vs_", "", model_key),
                    test_name = res$test_name,
                    statistic_name = stat_label,
                    statistic_value = res$statistic_value,
                    degrees_of_freedom = df_report,
                    n_samples = res$n_samples,
                    p_value = res$p_value,
                    adjusted_p_value = res$adjusted_p_value,
                    significance = res$significance,
                    is_significant = res$is_significant,
                    effect_size_r = res$effect_size_r,
                    conf_level = res$conf_level,
                    conf_int_low = res$conf_int_low,
                    conf_int_high = res$conf_int_high,
                    apa_report = apa,
                    apa_report_strict = apa_strict
                )
            }
        )

        return(p_res)
    }

    # Create the plot
    p <- ggplot(model_data, aes(x = model_display, y = r2, fill = model_class)) +
        # Add boxplots
        geom_boxplot(alpha = 0.7, width = 0.5) +    
        # Customize colors and labels
        scale_fill_manual(values = class_colors, name = "Model Class")

    if (show_connected_lines) {
        # Add lines connecting the same cell_type across models
        for (cell in unique(model_data$cell_type)) {
            for (org in unique(model_data$organism)) {
                cell_data <- model_data %>%
                    filter(cell_type == cell, organism == org)

                if (nrow(cell_data) > 1) {
                    p <- p + geom_line(
                        data = cell_data,
                        color = line_color,
                        linetype = linetype,
                        aes(group = interaction(cell_type, organism), color = organism),
                        alpha = line_alpha, linewidth = line_width
                    )
                }
            }
        }
        p <- p +
            # Add points for each observation
            geom_point(aes(color = organism), size = point_size, alpha = 0.8) +
            scale_color_manual(values = organism_colors, name = "Organism")
    } else {
        p <- p +
            # Add points for each observation
            ggforce::geom_sina(aes(color = organism), size = point_size, alpha = 0.8) +
            scale_color_manual(values = organism_colors, name = "Organism")
    }

    # Add cell type labels if requested
    if (show_cell_types) {
        p <- p + ggrepel::geom_text_repel(
            aes(label = cell_type),
            size = 3,
            box.padding = 0.5,
            point.padding = 0.3,
            segment.color = "grey50",
            min.segment.length = 0.1,
            max.overlaps = 15
        )
    }

    if (show_signif) {
        # Add significance bars and annotations
        for (i in 1:length(pairwise_results)) {
            result <- pairwise_results[[i]]

            # Only add if significant or we're not hiding non-significant results
            if (!hide_ns_bars || result$is_significant) {
                p <- p +
                    # Add significance line
                    annotate("segment",
                        x = result$position_a,
                        xend = result$position_b,
                        y = result$annotation_y,
                        yend = result$annotation_y
                    )

                if (result$is_significant) {
                    if (show_exact_pval) {
                        label <- paste0("p = ", format(result$adjusted_p_value, digits = 3), " ", result$significance)
                    } else {
                        label <- result$significance
                    }
                    p <- p +
                        # Add p-value text
                        annotate("text",
                            x = (result$position_a + result$position_b) / 2,
                            y = result$annotation_y * 1.02,
                            label = label,
                            size = 2
                        )
                } else {
                    p <- p +
                        # Add p-value text
                        annotate("text",
                            x = (result$position_a + result$position_b) / 2,
                            y = result$annotation_y * 1.02,
                            label = "ns",
                            size = 3
                        )
                }
            }
        }
    }
        
    # Handle subtitle
    if (is.null(subtitle)) {
        subtitle <- paste0(test_name, " with ", p_adjust_method, " correction")
        if (!is.null(explicit_comparisons)) {
            subtitle <- paste0(subtitle, " (explicit comparisons)")
        } else if (test_all_pairs) {
            subtitle <- paste0(subtitle, " (all pairs tested)")
        }
    }

    # Complete the plot
    p <- p +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +        
        # Labels
        labs(
            title = title,
            subtitle = subtitle,
            x = "",
            y = "R²"
        ) +
        ylim(ylim_range)

    # Return the plot object
    return(p)
}

# Function to create scatter plots comparing observed vs model predictions
plot_obs_model_scatter <- function(traj_preds, x_axis = "obs", x_axis_alias = x_axis, y_axis = "IQ", y_axis_alias = y_axis, xlab = x_axis_alias, ylab = y_axis_alias, point_size = 0.5, plot_r2 = TRUE, norm = TRUE){
    traj_preds <- as_tibble(traj_preds)
    
    # Create aliases for x and y axes if provided
    if (!is.null(x_axis_alias)) {
        traj_preds[[x_axis_alias]] <- traj_preds[[x_axis]]
        x_axis <- x_axis_alias
    }
    
    if (!is.null(y_axis_alias)) {
        traj_preds[[y_axis_alias]] <- traj_preds[[y_axis]]
        y_axis <- y_axis_alias
    }

    test_preds <- traj_preds %>%
        filter(!tss, type == "test")
    
    if (norm) {
        test_preds[[y_axis_alias]] <- norm01(test_preds[[y_axis_alias]])
        test_preds[[x_axis_alias]] <- norm01(test_preds[[x_axis_alias]])
        test_preds[[y_axis_alias]] <- rescale(test_preds[[y_axis_alias]], test_preds$obs)
        test_preds[[x_axis_alias]] <- rescale(test_preds[[x_axis_alias]], test_preds$obs)
    }
    
    r2 <- round(cor(test_preds[[x_axis]], test_preds[[y_axis]])^2, 3)
    t2_text <- paste0("R² = ", r2)

    p <- test_preds %>%
        ggplot(aes(x = !!sym(x_axis), y = !!sym(y_axis))) +
        geom_dense_scatter(size = point_size) +
        labs(x = xlab, y = ylab) +
        geom_abline() +
        xlim(-1, 1) +
        ylim(-1, 1) +
        theme(aspect.ratio = 1) +
        theme(legend.position = "none")
    if (plot_r2) {
        p <- p +
            labs(subtitle = t2_text)
    }
    return(p)
}
