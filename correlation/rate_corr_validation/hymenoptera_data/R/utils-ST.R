
get_branch_rate_across_smaps <- function(smaps) {
  total_changes <- 0
  branch.length <- smaps[[1]]$edge.length
  # i=2
  for (i in seq_along(smaps)) {
    S <- smaps[[i]]
    n.changes <- lapply(S$maps, function(x) length(x) - 1) %>% unlist()
    total_changes <- total_changes + n.changes
  }
  total_changes <- total_changes/length(smaps)
  total_changes <- total_changes/branch.length
  return(total_changes)
}

get_branch_rate_sample_smaps <- function(smaps) {
  total_changes <- vector('list', length = length(smaps))
  branch.length <- smaps[[1]]$edge.length
  # i=1
  for (i in seq_along(smaps)) {
    S <- smaps[[i]]
    n.changes <- lapply(S$maps, function(x) length(x) - 1) %>% unlist()
    total_changes[[i]] <- n.changes/branch.length
  }
  return(total_changes)
}

trim_xy_quantiles <- function(x, y, lower = 0.05, upper = 0.95, remove_zeros = FALSE, scale = FALSE, log.trans=FALSE, log.constant=0) {
  if (length(x) != length(y)) stop("Vectors x and y must be of equal length.")
  
  # Combine into data.frame for easy filtering
  df <- data.frame(x = x, y = y)
  
  # Optionally remove zeros
  if (remove_zeros) {
    df <- df[df$x != 0 & df$y != 0, ]
  }
  
  # Define quantile thresholds
  xq <- quantile(df$x, probs = c(lower, upper), na.rm = TRUE)
  yq <- quantile(df$y, probs = c(lower, upper), na.rm = TRUE)
  
  # Keep only values within both quantile ranges
  keep <- with(df, x >= xq[1] & x <= xq[2] & y >= yq[1] & y <= yq[2])
  df <- df[keep, ]
  
  # Optionally apply Z-transform
  if (scale) {
    df$x <- scale(df$x)[, 1]
    df$y <- scale(df$y)[, 1]
  }
  
  # Log transform
  if (log.trans) {
    df$x <- log(df$x + log.constant)
    df$y <- log(df$y + log.constant)
  }
  
  return(list(x = df$x, y = df$y))
}


test_correlation_with_resampling <- function(
    smap1, smap2,
    n.samples = 1000,
    cor_method = "spearman",
    remove_zeros = FALSE,
    scale = FALSE,
    log.trans = TRUE,
    log.constant = 1e-6,
    quantile.trim = c(0.05, 0.95)
) {
  # --- Observed correlation ---
  x_obs <- get_branch_rate_across_smaps(smap1)
  y_obs <- get_branch_rate_across_smaps(smap2)
  
  xy_obs <- trim_xy_quantiles(x_obs, y_obs,
                              lower = quantile.trim[1], upper = quantile.trim[2],
                              remove_zeros = remove_zeros, scale = scale,
                              log.trans = log.trans, log.constant = log.constant)
  
  obs_corr <- cor(xy_obs$x, xy_obs$y, method = cor_method, use = "pairwise.complete.obs")
  
  # --- Null distribution ---
  X <- get_branch_rate_sample_smaps(smap1)
  Y <- get_branch_rate_sample_smaps(smap2)
  
  n.maps <- length(X)
  n.branches <- length(X[[1]])
  
  null.corr <- numeric(n.samples)
  
  for (i in seq_len(n.samples)) {
    # Resample maps and branches
    X_resampled <- lapply(sample(X, size = n.maps, replace = TRUE), function(x) sample(x, n.branches, replace = TRUE))
    Y_resampled <- lapply(sample(Y, size = n.maps, replace = TRUE), function(x) sample(x, n.branches, replace = TRUE))
    
    x <- Reduce("+", X_resampled) / n.maps
    y <- Reduce("+", Y_resampled) / n.maps
    
    xy <- trim_xy_quantiles(x, y,
                            lower = quantile.trim[1], upper = quantile.trim[2],
                            remove_zeros = remove_zeros, scale = scale,
                            log.trans = FALSE, log.constant = 0)
    
    null.corr[i] <- cor(xy$x, xy$y, method = cor_method, use = "pairwise.complete.obs")
  }
  
  # --- P-value: two-tailed
  pval <- mean(abs(null.corr) >= abs(obs_corr), na.rm = TRUE)
  
  return(list(
    obs_corr = obs_corr,
    pval = pval,
    null_distribution = null.corr
  ))
}


#stm_list <- stm_merg
run_all_pairwise_correlations <- function(
    stm_list,
    n.samples = 1000,
    cor_method = "spearman",
    remove_zeros = FALSE,
    scale = FALSE,
    log.trans = TRUE,
    log.constant = 1e-6,
    quantile.trim = c(0.05, 0.95)
) {
  # Get all unique pairwise combinations
  parts <- names(stm_list)
  combos <- combn(parts, 2, simplify = FALSE)
  
  # Store results in a named list
  results <- list()
  
  for (pair in combos) {
    part1 <- pair[1]
    part2 <- pair[2]
    cat("Running correlation:", part1, "vs", part2, "\n")
    
    res <- tryCatch({
      test_correlation_with_resampling(
        stm_list[[part1]],
        stm_list[[part2]],
        n.samples = n.samples,
        cor_method = cor_method,
        remove_zeros = remove_zeros,
        scale = scale,
        log.trans = log.trans,
        log.constant = log.constant,
        quantile.trim = quantile.trim
      )
    }, error = function(e) {
      warning(sprintf("Failed for %s vs %s: %s", part1, part2, e$message))
      return(NULL)
    })
    
    res[['pair']] <- pair
    results[[paste(part1, part2, sep = " vs ")]] <- res
  }
  
  return(results)
}

run_all_pairwise_correlations_parallel <- function(
    stm_list,
    n.samples = 1000,
    cor_method = "spearman",
    remove_zeros = FALSE,
    scale = FALSE,
    log.trans = TRUE,
    log.constant = 1e-6,
    quantile.trim = c(0.05, 0.95),
    n.cores = parallel::detectCores() - 1
) {
  # Get all unique pairwise combinations
  parts <- names(stm_list)
  combos <- combn(parts, 2, simplify = FALSE)
  
  # Function to run correlation for one pair
  run_pair <- function(pair) {
    part1 <- pair[1]
    part2 <- pair[2]
    message("Running: ", part1, " vs ", part2)
    
    res <- tryCatch({
      test_correlation_with_resampling(
        stm_list[[part1]],
        stm_list[[part2]],
        n.samples = n.samples,
        cor_method = cor_method,
        remove_zeros = remove_zeros,
        scale = scale,
        log.trans = log.trans,
        log.constant = log.constant,
        quantile.trim = quantile.trim
      )
    }, error = function(e) {
      warning(sprintf("Failed for %s vs %s: %s", part1, part2, e$message))
      return(NULL)
    })
    
    res$pair <- pair
    return(res)
  }
  
  # Parallel execution using mclapply (UNIX/macOS)
  results_list <- parallel::mclapply(combos, run_pair, mc.cores = n.cores)
  
  # Name results for clarity
  names(results_list) <- sapply(combos, function(x) paste(x, collapse = " vs "))
  
  return(results_list)
}  


pval <- function(list, use.bonf=FALSE, alpha=NULL){
  if (use.bonf){
    pval <- unlist(lapply(list, function(x) x$pval))
    pval <- p.adjust(pval, method = "bonferroni")
  } else {
    pval <- unlist(lapply(list, function(x) x$pval))
  }
  
  
  if (!is.null(alpha)){
    pval <- pval[pval<=alpha]
  }
  return(pval)
}



corr_val <- function(list){
  unlist(lapply(list, function(x) x$obs_corr))
}

entropy <- function(p, epsilon = 1e-6, normalize = TRUE) {
  # Replace zeros with a small constant
  p[p == 0] <- epsilon
  
  # Normalize the vector to sum to 1
  p <- p / sum(p)
  
  # Compute raw entropy using natural logarithm
  raw_entropy <- -sum(p * log(p))
  
  if (normalize) {
    max_entropy <- log(length(p))  # Max entropy with same number of categories
    return(raw_entropy / max_entropy)
  } else {
    return(raw_entropy)
  }
}

get_enrichment_mask <- function(rate_matrix, lower_q = 0.025, upper_q = 0.975,
                                direction = c("both", "over", "under"),
                                return_type = c("binary", "masked")) {
  
  direction <- match.arg(direction)
  return_type <- match.arg(return_type)
  
  enrichment_mask <- matrix(0, nrow = nrow(rate_matrix), ncol = ncol(rate_matrix))
  colnames(enrichment_mask) <- colnames(rate_matrix)
  rownames(enrichment_mask) <- rownames(rate_matrix)
  
  for (col in colnames(rate_matrix)) {
    col_vals <- rate_matrix[, col]
    q <- quantile(col_vals, probs = c(lower_q, upper_q), na.rm = TRUE)
    
    if (direction == "both") {
      mask <- col_vals < q[1] | col_vals > q[2]
    } else if (direction == "over") {
      mask <- col_vals > q[2]
    } else {  # under
      mask <- col_vals < q[1]
    }
    
    if (return_type == "binary") {
      enrichment_mask[, col] <- as.numeric(mask)
    } else {  # masked
      enrichment_mask[, col] <- ifelse(mask, col_vals, 0)
    }
  }
  
  return(enrichment_mask)
}


plot_branch_colored_tree <- function(tree, branch_values, palette = viridis::viridis(100), 
                                     title = "Branch Coloring", legend_title = "Value") {
  if (length(branch_values) != Nedge(tree)) {
    stop("Length of branch_values must match the number of branches in the tree.")
  }
  
  # Scale values to palette
  value_scaled <- as.numeric(cut(branch_values, breaks = length(palette), labels = FALSE))
  branch_colors <- palette[value_scaled]
  
  # Plot tree with colored branches
  plot(tree, edge.color = branch_colors, edge.width = 2, cex = 0.5)
  title(title)
  
  # Create legend labels
  min_val <- round(min(branch_values, na.rm = TRUE), 2)
  max_val <- round(max(branch_values, na.rm = TRUE), 2)
  mid_val <- round(mean(range(branch_values, na.rm = TRUE)), 2)
  
  # Legend colors
  low_col  <- palette[1]
  mid_col  <- palette[round(length(palette) / 2)]
  high_col <- palette[length(palette)]
  
  # Add legend
  legend("topright",
         legend = c(paste0("Low (", min_val, ")"),
                    paste0("Medium (", mid_val, ")"),
                    paste0("High (", max_val, ")")),
         fill = c(low_col, mid_col, high_col),
         border = NA,
         cex = 0.7,
         title = legend_title)
}

