library(phytools)
library(parallel)

count_followed <- function(vec, a, b) {
  # Convert to character if needed, but works for both numeric and character
  vec <- as.character(vec)
  a <- as.character(a)
  b <- as.character(b)
  
  # Handle short vectors
  if (length(vec) < 2) return(0)
  
  # Count occurrences where a is followed by b
  count <- sum(vec[-length(vec)] == a & vec[-1] == b)
  return(count)
}

count_focal_changes <- function(simmap, focal.changes) {
  # Helper to count how many times "from" → "to" occurs in a vector of states
  count_followed <- function(vec, from, to) {
    if (length(vec) < 2) return(0)
    sum(vec[-length(vec)] == from & vec[-1] == to)
  }
  
  # Parse focal changes into from→to pairs
  focal_pairs <- strsplit(focal.changes, ">")
  n_maps <- length(simmap)
  
  # Get number of branches from first simmap
  n_branches <- length(simmap[[1]]$maps)
  branch_names <- names(simmap[[1]]$maps)
  
  # Results array: maps × branches × focal changes
  results_array <- array(
    0,
    dim = c(n_maps, n_branches, length(focal_pairs)),
    dimnames = list(
      paste0("map", seq_len(n_maps)),
      branch_names,
      focal.changes
    )
  )
  
  # Loop over all simmap replicates
  for (m in seq_len(n_maps)) {
    map.changes <- lapply(simmap[[m]]$maps, names)
    
    # Loop over focal changes
    for (f in seq_along(focal_pairs)) {
      from <- focal_pairs[[f]][1]
      to   <- focal_pairs[[f]][2]
      
      results_array[m, , f] <- unlist(
        lapply(map.changes, function(x) count_followed(x, from, to))
      )
    }
  }
  
  # Per-branch totals across all maps
  per_branch <- apply(results_array, c(2, 3), sum)
  
  # Overall totals per focal change
  total_summary <- colSums(per_branch)
  
  return(list(
    per_map = results_array,     # detailed counts per map × branch × transition
    per_branch = per_branch,     # summed across all maps per branch
    total = total_summary        # final overall counts
  ))
}


get_per_branch_rate <- function(simmap, focal.changes){
  res <- count_focal_changes(simmap, focal.changes)
  apply(res$per_branch, 1, sum) / length(simmap) / simmap[[1]]$edge.length
}


library(phytools)
library(parallel)

simulate_multi_char_simmap <- function(fit_Q, tree, char.mt, dt, p.root, Nsimmap = 100, ncores = detectCores() - 1) {
  
  # Fix Q matrix (CTMC valid)
  Q <- fit_Q$solution
  Q_fixed <- Q
  Q_fixed[is.na(Q_fixed)] <- 0
  diag(Q_fixed) <- -rowSums(Q_fixed)
  
  # Get all character names
  char_names <- names(char.mt)
  
  # Start parallel cluster
  cl <- makeCluster(ncores)
  
  # Export needed objects to all workers
  clusterExport(cl, varlist = c("char.mt", "char_names", "tree", "Q_fixed", "p.root", "Nsimmap", "dt"), envir = environment())
  
  # Load required package on each worker
  clusterEvalQ(cl, library(phytools))
  
  # Run parallel simulation over characters
  simmap_multi_char <- parLapply(cl, char_names, function(char_name) {
    # Extract matrix for current character
    char <- char.mt[[char_name]] #[1:Ntip(tree), ]
    #rownames(char) <- dt[, 1]
    
    # Simulate stochastic maps
    make.simmap(tree, char, Q = Q_fixed, nsim = Nsimmap, pi = p.root)
  })
  
  # Stop cluster
  stopCluster(cl)
  
  # Name results
  names(simmap_multi_char) <- char_names
  
  return(simmap_multi_char)
}


sum_per_branch_rates <- function(simmap_results, direc.changes) {
  # Apply get_per_branch_rate to each element of simmap_results
  branch_rate_list <- lapply(simmap_results, function(x) get_per_branch_rate(x, direc.changes))
  
  # Ensure all elements are numeric vectors of the same length
  branch_rate_list <- lapply(branch_rate_list, function(x) as.numeric(x))
  
  # Sum all vectors element-wise
  branch_rate_sum <- Reduce(`+`, branch_rate_list)
  
  # return(list(
  #   per_char = branch_rate_list,  # individual vectors for each char
  #   total = branch_rate_sum       # summed vector across all chars
  # ))
  return(branch_rate_sum)
}


recode_simmap_states <- function(simmap, recode) {
  # Parse the recode rules into a named vector
  recode_pairs <- strsplit(recode, ">")
  from <- sapply(recode_pairs, `[`, 1)
  to   <- sapply(recode_pairs, `[`, 2)
  recode_map <- setNames(to, from)
  
  # Function to recode a single branch map
  recode_branch <- function(branch_map) {
    names(branch_map) <- ifelse(
      names(branch_map) %in% names(recode_map),
      recode_map[names(branch_map)],
      names(branch_map)
    )
    return(branch_map)
  }
  
  # Apply to each edge of each tree in the simmap object
  if (inherits(simmap, "multiSimmap")) {
    for (k in seq_along(simmap)) {
      simmap[[k]]$maps <- lapply(simmap[[k]]$maps, recode_branch)
    }
  } else if (inherits(simmap, "simmap")) {
    simmap$maps <- lapply(simmap$maps, recode_branch)
  } else {
    stop("Input must be a simmap or multiSimmap object.")
  }
  
  return(simmap)
}

