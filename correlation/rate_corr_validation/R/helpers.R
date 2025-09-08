
# Simulate informative tip data.

sim_info_data <- function(tree, focal_edge, Q_base, Q_jump, nsim) {
  
  direct_data <- c()

    for (k in 1:nsim) {
    
    # Loop until simulate an informative trait.
    x <- c()
    while (length(unique(x)) != 2) {
      
      x <- simDirectional(tree, focal_edge, Q_base, Q_jump, nsim = 1)
      
    }
    
    direct_data <- cbind(direct_data, x)
    
    }
  
  colnames(direct_data) <- paste0("C", 1:nsim)
  
  return(direct_data)
  
}


# Sample stochastic histories conditioned on data.

sample_cond <- function(tree, tip_data, nsim, nstm) {
  
  # List of fitted models.
  fitted_models <- setNames(vector(mode = "list", length = nsim), paste0("C", 1:nsim))
  
  # List of stochastic maps.
  smaps <- setNames(vector(mode = "list", length = nsim), paste0("C", 1:nsim))
  
  i = 1
  # Loop over all traits.
  for (i in 1:nsim){
    
    cat(paste0("\n", "Working on: ", "C", i, ": ", Sys.time(), "\n"))
    
    # Set character vector.
    char <- cbind(rownames(tip_data), tip_data[,i])
    
    # Build set of models.
    mods <- c("ER", "ARD")
    
    # Fit models.
    fit_Q <- setNames(vector(mode = "list", length = length(mods)), mods)
    
    j = 1
    for (j in 1:length(mods)) {
      
      fit_Q[[j]] <- corHMM(phy = tree, data = char, model = mods[[j]], rate.cat = 1, root.p = "maddfitz")
      
    }
    
    # Get best model.
    w <- geiger::aicw(sapply(fit_Q, function(x) x$AICc))[,3]
    
    # Get Q matrix.
    Q <- fit_Q[[min(which(w == max(w)))]]$solution
    
    # Store fitted Q matrix.
    fitted_models[[i]] <- Q
    
    # Sample stochastic maps.
    smaps[[i]] <- makeSimmap(tree = tree, data = char, model = Q, rate.cat = 1, nSim = nstm)
    
  }
  
  # Return results.
  return(list(smaps = smaps, models = fitted_models))
  
}


# Sample stochastic histories unconditioned on data.

sim_uncond <- function(tree, fitted_models, nsim, nstm) {
  
  # List of stochastic maps.
  smaps <- setNames(vector(mode = "list", length = nsim), paste0("C",1:nsim))
  
  for (i in 1:nsim) {
    
    cat(paste0("\n", "Working on: ", "C", i, ": ", Sys.time(), "\n"))
    
    # Build Q matrix from models fitted to the observed data.
    Q <- fitted_models[[i]]
    diag(Q) <- -rowSums(Q, na.rm = T)
    Q[is.na(Q)] <- 0
    colnames(Q) <- rownames(Q) <- 0:(dim(Q)[1] - 1)
    
    # Simulate character histories.
    smaps[[i]] <- sim.history(tree, Q = Q, nsim = nstm)
    
  }
  
  # Return results.
  return(smaps)
  
}

# Calculate total Hamming distance per branch.

hm_dist_br <- function(stm) {
  
  br_hm_dist <- vector(mode = "list", length = length(stm))
  
  for (i in 1:length(stm)) {
    
    br_hm <- sapply(stm[[i]]$maps, function(x) sum(stringdist::stringdist(tail(names(x), -1), head(names(x), -1), method = "hamming")) )
    br_hm[is.infinite(br_hm)] <- 0
    br_hm_dist[[i]] <- br_hm
    
  }
  
  return(br_hm_dist)
  
} 

# Permutation test.

# x <- hm_dist_m[,"A1"]
# y <- hm_dist_m[,"B1"] 
# nperm = 10
# method = "spearman"

perm_test <- function(x, y, nperm = 10000, method = c("pearson", "spearman")) {
  
  # Get the observed correlation.
  obs_corr <- cor(x, y, method = method)
  
  # Set a vector to store correlation values.
  perm_corr <- numeric(nperm)
  
  # Permutations.
  # i=1
  for (i in 1:nperm) { perm_corr[[i]] <- cor(x, sample(y, size = length(y), replace = F), method = method) }
  
  # Get p-values.
  if (obs_corr > 0) { 
    
    pval <- (sum(perm_corr >= obs_corr)/nperm)
    
    return(list("obs_corr" = obs_corr, "pval" = pval, "perm_corr" = perm_corr) ) 
    
  } else { 
    
    pval = (sum(perm_corr <= obs_corr)/nperm)
    
    return(return(list("obs_corr" = obs_corr, "pval" = pval, "perm_corr" = perm_corr) ))
    
  }
  
}



perm_test_Double <- function(x, y, nperm = 10000, method = c("pearson", "spearman")) {
  
  # Get the observed correlation.
  obs_corr <- cor(x, y, method = method)
  
  # Set a vector to store correlation values.
  perm_corr <- numeric(nperm)
  
  # Permutations.
  # i=1
  for (i in 1:nperm) {
    y.perm <- sample(y, size = length(y), replace = F)
    x.perm <- sample(x, size = length(x), replace = F)
    perm_corr[[i]] <- cor(x, sample(y, size = length(y), replace = F), method = method) 
    }
  
  # Get p-values.
  if (obs_corr > 0) { 
    
    pval <- (sum(perm_corr >= obs_corr)/nperm)
    
    return(list("obs_corr" = obs_corr, "pval" = pval, "perm_corr" = perm_corr) ) 
    
  } else { 
    
    pval = (sum(perm_corr <= obs_corr)/nperm)
    
    return(return(list("obs_corr" = obs_corr, "pval" = pval, "perm_corr" = perm_corr) ))
    
  }
  
}