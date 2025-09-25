
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

# This function test correlations and also perfoms
# - quantile.trim: outlier removal by quantile trimming
# - remove_zeros: removing data where at least one branch is zero
# - log.trans: log transform (log.constant is added to zero values as log(0) is INF)
# - scale: Z transform
test_correlation_with_resampling <- function(
    smap1, smap2,
    n.samples = 1000,
    cor_method = "spearman",
    remove_zeros = FALSE,
    scale = FALSE,
    log.trans = TRUE,
    log.constant = 1e-6,
    quantile.trim = c(0.05, 0.95),
    test=NULL
) {
  # --- Observed correlation ---
  if (is.null(test)){
    x_obs <- get_branch_rate_across_smaps(smap1)
    y_obs <- get_branch_rate_across_smaps(smap2)
    # --- Null distribution ---
    X <- get_branch_rate_sample_smaps(smap1)
    Y <- get_branch_rate_sample_smaps(smap2)
  } else {
    # x_obs <- test$x_obs
    # y_obs <- test$y_obs
    x_obs <- Reduce("+", test$x_obs) / length(test$x_obs)
    y_obs <- Reduce("+", test$y_obs) / length(test$y_obs)
    X <-test$x_obs
    Y <- test$y_obs
  } 

  
  xy_obs <- trim_xy_quantiles(x_obs, y_obs,
                              lower = quantile.trim[1], upper = quantile.trim[2],
                              remove_zeros = remove_zeros, scale = scale,
                              log.trans = log.trans, log.constant = log.constant)
  
  obs_corr <- cor(xy_obs$x, xy_obs$y, method = cor_method, use = "pairwise.complete.obs")
  
  # --- Null distribution ---
  # X <- get_branch_rate_sample_smaps(smap1)
  # Y <- get_branch_rate_sample_smaps(smap2)
  
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
    xy_obs = xy_obs,
    null_distribution = null.corr
  ))
}


test_correlation_R <- function(
    smap1, smap2,
    cor_method = "spearman",
    remove_zeros = FALSE,
    scale = FALSE,
    log.trans = TRUE,
    log.constant = 1e-6,
    quantile.trim = c(0.05, 0.95),
    test=NULL
) {
  # --- Observed correlation ---
  if (is.null(test)){
    x_obs <- get_branch_rate_across_smaps(smap1)
    y_obs <- get_branch_rate_across_smaps(smap2)
  } else {
    x_obs <- Reduce("+", test$x_obs) / length(test$x_obs)
    y_obs <- Reduce("+", test$y_obs) / length(test$y_obs)
  } 
  
  
  xy_obs <- trim_xy_quantiles(x_obs, y_obs,
                              lower = quantile.trim[1], upper = quantile.trim[2],
                              remove_zeros = remove_zeros, scale = scale,
                              log.trans = log.trans, log.constant = log.constant)
  
  #obs_corr <- cor(xy_obs$x, xy_obs$y, method = cor_method, use = "pairwise.complete.obs")
  obs_corr <- cor.test(xy_obs$x, xy_obs$y, alternative = "two.sided", method = cor_method, use = "pairwise.complete.obs")
  #obs_corr <- cor.test(x, y, method = "pearson", use = "pairwise.complete.obs", alternative = "two.sided")
  
  return(list(
    obs_corr = obs_corr$estimate,
    pval = obs_corr$p.value,
    xy_obs = xy_obs
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
    use.resampling=TRUE,
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
      if (use.resampling){
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
      } else {
        test_correlation_R(
          stm_list[[part1]],
          stm_list[[part2]],
          cor_method = cor_method,
          remove_zeros = remove_zeros,
          scale = scale,
          log.trans = log.trans,
          log.constant = log.constant,
          quantile.trim = quantile.trim
        )
      }
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



# corr_val <- function(list){
#   unlist(lapply(list, function(x) x$obs_corr))
# }

corr_val <- function(list){
  vals <- sapply(list, function(x) x$obs_corr)
  names(vals) <- names(list)
  return(vals)
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



create_pairwise_matrix <- function(cv, pv, dist_matrix_dist, mode = c("na", "zero", "remove", 'none')) {
  mode <- match.arg(mode)

  # Extract unique part names
  parts <- unique(unlist(strsplit(names(cv), " vs ")))
  parts <- sort(parts)  # ensure consistent ordering

  # Initialize correlation matrix
  cor_mat <- matrix(NA, nrow = length(parts), ncol = length(parts),
                    dimnames = list(parts, parts))

  # Fill in correlation values
  for (pair in names(cv)) {
    split_pair <- unlist(strsplit(pair, " vs "))
    i <- split_pair[1]
    j <- split_pair[2]

    if (pv[pair] < 0.05) {
      cor_mat[i, j] <- abs(cv[pair])
      cor_mat[j, i] <- abs(cv[pair])
    } else {
      if (mode == "zero" || mode == "remove") {
        cor_mat[i, j] <- 0
        cor_mat[j, i] <- 0
      } else if (mode == "na") {
        cor_mat[i, j] <- NA
        cor_mat[j, i] <- NA
      }
    }
  }

  # Convert dist_matrix_dist to full square matrix and match names
  dist_mat <- as.matrix(dist_matrix_dist)
  all_parts <- intersect(rownames(cor_mat), rownames(dist_mat))
  dist_mat <- dist_mat[all_parts, all_parts]
  cor_mat <- cor_mat[all_parts, all_parts]

  # Handle removal of rows/cols for 'remove' mode
  if (mode == "remove") {
    keep <- apply(cor_mat, 1, function(row) any(row != 0 & !is.na(row)))
    cor_mat <- cor_mat[keep, keep]
    dist_mat <- dist_mat[keep, keep]
  }

  # Return distance objects + matrices with matched order
  cor_mat_dist <- as.dist(cor_mat)
  dist_matrix_dist_out <- as.dist(dist_mat)

  return(list(
    #cor_mat = cor_mat,
    #dist_matrix = dist_mat,
    cor_mat_dist = cor_mat_dist,
    dist_mat_dist = dist_matrix_dist_out
  ))
}



# Suppose you already have:
# dist_matrix <- distances(g)
# ana.d <- c("cranium vs pronotum", "cranium vs propectus", ...)

# Function to extract distances
extract_pairs <- function(dist_matrix, pairs) {
  results <- numeric(length(pairs))
  
  for (i in seq_along(pairs)) {
    parts <- strsplit(pairs[i], " vs ")[[1]]
    results[i] <- dist_matrix[parts[1], parts[2]]
  }
  
  names(results) <- pairs
  return(results)
}


get_len <- function(list_elem){
  unlist(lapply(list_elem, function(x) length(x) ))
}



plot_enrichment_stacked <- function(tree, enr, n_points = 100, denom = c("active", "enriched")) {
  # enr: matrix (branches × body regions), 0/1 enrichment
  denom <- match.arg(denom)  # choose between "active" or "enriched"
  
  # branch times
  H <- nodeHeights(tree)
  max_height <- max(H)
  times <- seq(0, max_height, length.out = n_points)
  
  regions <- colnames(enr)
  props <- matrix(0, nrow = n_points, ncol = length(regions))
  colnames(props) <- regions
  
  #i=10
  for (i in seq_along(times)) {
    t <- times[i]
    edge_ids <- which(t >= H[,1] & t <= H[,2])
    if (length(edge_ids) > 0) {
      enr_now <- enr[edge_ids, , drop = FALSE]
      if (denom=='active'){
        props[i, ] <- colSums(enr_now) / length(edge_ids)
      } else {
        props[i, ] <- colSums(enr_now) / ncol(enr) 
      }
    }
  }
  
  # make into data.frame
  df <- data.frame(time = max_height - times, props, check.names = FALSE)
  
  # melt to long format
  df_long <- reshape2::melt(df, id.vars = "time",
                            variable.name = "region",
                            value.name = "prop")
  
  ggplot(df_long, aes(x = time, y = prop, fill = region)) +
    geom_area(position = "stack") +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_x_reverse(expand = c(0, 0)) +  # root left, present right
    labs(x = "Time before present", y = "Proportion of enriched lineages") +
    theme_classic()
}




# states <- c("10", "00", "00")
# states <- c("00", '01')
# states <- c('01', '01')
# diversity_metric(states, type = "weighted")
# diversity_metric(states, type = "Neff_scaled")

diversity_metric <- function(states, type = c("Neff", "Neff_scaled", "weighted"), base = exp(1)) {
  type <- match.arg(type)
  
  shannon_entropy <- function(x) {
    freqs <- table(x) / length(x)
    -sum(freqs * log(freqs, base = base))
  }
  
  if (type == "Neff") {
    enriched_states <- states[grepl("1", states)]
    if (length(enriched_states) == 0) return(0)
    H_plus <- shannon_entropy(enriched_states)
    return(base ^ H_plus)
  }
  
  if (type == "Neff_scaled") {
    enriched_states <- states[grepl("1", states)]
    prop_enriched <- mean(grepl("1", states))
    if (length(enriched_states) == 0) return(0)
    H_plus <- shannon_entropy(enriched_states)
    Neff_plus <- base ^ H_plus
    return(Neff_plus * prop_enriched)
  }
  
  if (type == "weighted") {
    H <- shannon_entropy(states)
    f_enriched <- mean(grepl("1", states))
    return(H + f_enriched)
  }
}


#enr <- enr[,c(1:2)]
# enr <- enr[,1, drop=F]
#n_points = 10
plot_enrichment_entropy <- function(tree, enr, n_points = 100, time_points=NULL,  type=c("Neff", "Neff_scaled", "weighted"), base = exp(1)) {
  # enr: matrix (branches × body regions), 0/1 enrichment
  
  # branch times (start, end for each edge)
  H <- nodeHeights(tree)
  max_height <- max(H)
  if (is.null(time_points)){
    times <- seq(0, max_height, length.out = n_points)
    times[n_points] <- times[n_points] - 1e-4
  } else { # time points provided by user
    times <- time_points
  }

  
  #entropy_vals <- numeric(n_points)
  entropy_vals <- numeric(length(times))
  # i=10
  for (i in seq_along(times)) {
    t <- times[i]
    edge_ids <- which(t >= H[,1] & t <= H[,2])  # edges active at time t
    if (length(edge_ids) > 0) {
      # enrichment states for active branches
      enr_now <- enr[edge_ids, , drop = FALSE]
      states <- apply(enr_now, 1, paste0, collapse = "")
      # states
      entropy_vals[i] <- diversity_metric(states, type=type, base=base) 
      
    } else {
      entropy_vals[i] <- NA
    }
  }
  
  df <- data.frame(time = rev(times), entropy = entropy_vals)
  
  p <- ggplot(df, aes(x = time, y = entropy)) +
    geom_line() +
    theme_minimal() +
    scale_x_reverse(breaks = seq(0, max_height, by = 50)) +
    labs(x = "Time", y = "Entropy of enrichment")

  return(list(data = df, plot = p))
}



# states <- c("10", "00", "00", "01")
# states <- c("00", '00')
# relative_null_metric(c("00", '01'))
# relative_null_metric(states)

relative_null_metric <- function(states, type = c("not_scaled","scaled")) {
  type <- match.arg(type)
  
  if (length(states) == 0) return(NA)
  
  len <- nchar(states[1])              # number of body regions
  enriched <- states[grepl("1", states)]
  n <- length(enriched)
  if (n == 0) return(0)                # no enriched states
  
  # Convert enriched states to binary matrix
  mat <- strsplit(enriched, split = "")
  mat <- do.call(rbind, mat)
  mat <- apply(mat, 2, as.integer)
  if (is.null(dim(mat))) {             # only one enriched state
    mat <- matrix(mat, nrow = 1)
  }
  
  # Distances from null (00...0)
  d_null <- rowSums(mat)               # Hamming distance to all-zeros
  rel_dist <- mean(d_null/len)       # normalized mean distance
  
  if (type == "not_scaled") {
    return(rel_dist)
  }
  
  if (type == "scaled") {
    prop_enriched <- mean(grepl("1", states))
    return(rel_dist * prop_enriched)
  }
}


#type = c("scaled")
plot_enrichment_strength <- function(tree, enr, n_points = 100,  type = c("not_scaled","scaled")) {
  # enr: matrix (branches × body regions), 0/1 enrichment
  
  # branch times (start, end for each edge)
  H <- nodeHeights(tree)
  max_height <- max(H)
  times <- seq(0, max_height, length.out = n_points)
  times[n_points] <- times[n_points] - 1e-4
  
  entropy_vals <- numeric(n_points)
  # i=3
  for (i in seq_along(times)) {
    t <- times[i]
    edge_ids <- which(t >= H[,1] & t <= H[,2])  # edges active at time t
    if (length(edge_ids) > 0) {
      # enrichment states for active branches
      enr_now <- enr[edge_ids, , drop = FALSE]
      states <- apply(enr_now, 1, paste0, collapse = "")
      # states
      entropy_vals[i] <- relative_null_metric(states, type=type)
      
    } else {
      entropy_vals[i] <- NA
    }
  }
  
  df <- data.frame(time = rev(times), strength = entropy_vals)
  
  p <- ggplot(df, aes(x = time, y = strength)) +
    geom_line() +
    theme_minimal() +
    scale_x_reverse(breaks = seq(0, max_height, by = 50)) +
    labs(x = "Time", y = "Enrichment strength")
  
  return(list(data = df, plot = p))
}


enrich_next_level <- function(enr, region_class, group_order = NULL) {
  # enr: matrix branches × regions (0/1 enrichment)
  # region_class: data.frame with columns: region, group
  # group_order: optional vector giving desired order of groups
  
  # Convert to data.frame for easier handling
  enr_df <- as.data.frame(enr)
  enr_df$branch <- seq_len(nrow(enr_df))  # keep branch ID
  
  # Reshape to long format
  library(reshape2)
  enr_long <- melt(enr_df, id.vars = "branch",
                   variable.name = "region",
                   value.name = "enriched")
  
  # Join with region_class
  enr_long <- merge(enr_long, region_class, by = "region")
  
  # Collapse: for each branch × group, mark 1 if any region in group enriched
  enr_collapsed <- aggregate(enriched ~ branch + group, data = enr_long,
                             FUN = function(x) as.integer(any(x > 0)))
  
  # Reshape back to wide matrix: branches × groups
  enr_class <- reshape(enr_collapsed,
                       idvar = "branch",
                       timevar = "group",
                       direction = "wide")
  
  # Clean column names
  colnames(enr_class) <- sub("enriched\\.", "", colnames(enr_class))
  
  # Remove branch col, set rownames
  rownames(enr_class) <- enr_class$branch
  enr_class <- enr_class[, -1, drop = FALSE]
  
  # Reorder groups if order provided
  if (!is.null(group_order)) {
    enr_class <- enr_class[, group_order, drop = FALSE]
  }
  
  return(enr_class)
}


get_tip_path <- function(tree, tip) {
  # tip can be label or index
  if (is.character(tip)) {
    tip <- which(tree$tip.label == tip)
  }
  
  # all ancestors of the tip (including root)
  anc <- phangorn::Ancestors(tree, tip, type = "all")
  
  # edge indices where the child is the tip or one of its ancestors
  path_ids <- which(tree$edge[,2] %in% c(tip, anc))
  
  return(path_ids)
}


get_tip_enrichment <- function(tree, enr.masked){
  tree.height <- max(node.depth.edgelength(tree))
  ntip <- Ntip(tree)
  
  enr.tip <- matrix(NA, nrow=ntip, ncol=ncol(enr.masked))
  # dim(enr.tip)
  colnames(enr.tip) <- colnames(enr.masked)
  #br_index <- 1
  for (br_index in 1:ncol(enr.masked)){
    br <- enr.masked[,br_index]
    #tip_total_rate <- rep(NA,  ntip)
    #tip_index <- 1
    for (tip_index in 1:ntip){
      #tip_index <- 1
      #tree$tip.label[77]
      tip_path <- get_tip_path(tree, tip_index)
      tip_path_rates <- br[tip_path]
      # we scale it by tot time over path
      enr.tip[tip_index, br_index] <- sum(tip_path_rates)/tree.height
    }
  }
  return(enr.tip)
}


# Function: edge times for each branch
# Function: edge times from tips
# start = distance from tip where the branch starts
# end   = distance from tip where the branch ends
# tr <- rtree(5)   # random phylogeny with 5 tips
# edge_times(tr)
# plot(tr)
# axisPhylo()
# edgelabels()
edge_times_from_tips <- function(tree) {
  if (is.null(tree$edge.length)) {
    stop("Tree must have edge lengths (branch lengths).")
  }
  
  # node depths = distance from root
  node_depths <- node.depth.edgelength(tree)
  tree_height <- max(node_depths)
  
  edges <- tree$edge
  # time from tip = tree height - distance from root
  start_from_tip <- tree_height - node_depths[edges[,2]]  # branch starts at child
  end_from_tip   <- tree_height - node_depths[edges[,1]]  # branch ends at parent
  
  data.frame(
    edge  = seq_len(nrow(edges)),
    start = start_from_tip,
    end   = end_from_tip
  )
}


# cor_mat <- cor_trimmed(
#   rates,
#   method = "pearson",
#   log_transform = TRUE,
#   log_const = 1e-6,
#   trim_quantiles = c(0.05, 0.95)
# )
cor_trimmed <- function(rates, 
                        method = "pearson", 
                        log_transform = FALSE, 
                        log_const = 1e-6,
                        trim_quantiles = NULL,
                        use = "pairwise.complete.obs") {
  # rates: matrix (branches x body regions)
  # method: correlation method ("pearson", "spearman", etc.)
  # log_transform: TRUE/FALSE
  # log_const: small constant added before log
  # trim_quantiles: e.g. c(0.05, 0.95) → remove values outside this range per BR
  # use: NA handling, same as cor()
  
  mat <- rates
  
  # log transform if requested
  if (log_transform) {
    mat <- log(mat + log_const)
  }

  
  # trim outliers per column if requested
  if (!is.null(trim_quantiles)) {
    for (j in seq_len(ncol(mat))) {
      col_vals <- mat[, j]
      qs <- quantile(col_vals, probs = trim_quantiles, na.rm = TRUE)
      keep <- col_vals >= qs[1] & col_vals <= qs[2]
      mat[!keep, j] <- NA  # set outliers to NA
    }
  }
  
  # # log transform if requested
  # if (log_transform) {
  #   mat <- log(mat + log_const)
  # }
  
  # compute correlation
  cor_mat <- cor(mat, use = use, method = method)
  
  return(cor_mat)
}


# cor_results <- rolling_rate_correlation(rates, edge_times,window_size = 30, step_size = 2)
# rates <- rates[, c(1:3)]
rolling_rate_correlation <- function(rates, edge_times, window_size = 20, step_size = 5,
                                     trim_quantiles = NULL, log_transform = FALSE, use.abs.value=TRUE) {
  # rates: matrix [edges x body_regions]
  # edge_times: data.frame with cols edge, start, end (from tips)
  # window_size: width of sliding window in time units
  # step_size: how much to slide at each step
  
  # Ensure matching edges
  if (nrow(rates) != nrow(edge_times)) {
    stop("rates rows must match edges in edge_times")
  }
  
  # Midpoint time for each edge
  edge_times$mid <- (edge_times$start + edge_times$end) / 2
  
  max_t <- max(edge_times$end)
  min_t <- min(edge_times$start)
  
  # Window centers
  centers <- seq(min_t, max_t, by = step_size)
  
  # Store results
  results <- list()
  
  # c=0
  for (c in centers) {
    # Window bounds
    lower <- c - window_size/2
    upper <- c + window_size/2
    
    # Select edges whose midpoint falls inside window
    in_window <- edge_times$mid >= lower & edge_times$mid <= upper
    sub_rates <- rates[in_window, , drop = FALSE]
    
    if (nrow(sub_rates) > 2) {
      #cor_mat <- cor(sub_rates, use = "pairwise.complete.obs", method='pearson')
      cor_mat <- cor_trimmed(
        sub_rates,
        method = "pearson",
        log_transform = log_transform,
        log_const = 1e-6,
        trim_quantiles = trim_quantiles, #c(0.01, 0.99),
        use = "pairwise.complete.obs"
      )
      # use absolute values
      if (use.abs.value){
        cor_mat <- abs(cor_mat)
      }
    } else {
      cor_mat <- matrix(NA, ncol = ncol(rates), nrow = ncol(rates))
      colnames(cor_mat) <- colnames(rates)
      rownames(cor_mat) <- colnames(rates)
    }
    
    results[[as.character(c)]] <- cor_mat
  }
  
  return(results)
}


count_branches_sliding <- function(edge_times, window_size = 20, step_size = 5, plot = TRUE) {
  # Determine overall time range
  t_min <- min(edge_times$end, na.rm = TRUE)
  t_max <- max(edge_times$start, na.rm = TRUE)
  
  # Define sliding windows
  centers <- seq(t_min, t_max, by = step_size)
  results <- data.frame(
    time = centers,
    n_branches = NA_integer_
  )
  
  # Count branches per window
  for (i in seq_along(centers)) {
    w_start <- centers[i] - window_size / 2
    w_end   <- centers[i] + window_size / 2
    
    overlap <- with(edge_times,
                    (start >= w_start & start <= w_end) | 
                      (end   >= w_start & end   <= w_end) |
                      (start <= w_start & end   >= w_end))
    
    results$n_branches[i] <- sum(overlap, na.rm = TRUE)
  }
  
  # Plot if requested
  if (plot) {
    library(ggplot2)
    p <- ggplot(results, aes(x = time, y = n_branches)) +
      geom_line(size = 1, color = "black") +
      scale_x_reverse(
        breaks = seq(0, 300, by = 25),
        expand = c(0, 0)) +
      theme_classic(base_size = 14) +
      labs(x = "Time (Ma)", y = "Number of branches")
    print(p)
  }
  
  return(results)
}


cor_list_to_df <- function(cor_results) {
  df <- do.call(rbind, lapply(names(cor_results), function(t) {
    mat <- cor_results[[t]]
    if (all(is.na(mat))) return(NULL)
    as.data.frame(as.table(mat)) %>%
      rename(region1 = Var1, region2 = Var2, cor = Freq) %>%
      mutate(time = as.numeric(t))
  }))
  return(df)
}
