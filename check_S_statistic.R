
# Load packages.
library(ape)
library(phytools)
library(tidyverse)

# Import data.
stm_amalg <- readRDS("stm_amalg.RDS")
stm_amalg_null <- readRDS("stm_amalg_null.RDS")

# Merge state categories across branches.
stm_merg <- merge_tree_cat_list(stm_amalg)
stm_merg_null <- merge_tree_cat_list(stm_amalg_null)

# Get max Hamming distance per branch.
# Data.
br_hm_dist <- vector(mode = "list", length = length(stm_merg))
for (i in 1:length(stm_merg)) {
  
  br_hm <- sapply(stm_merg[[i]]$maps, function(x) stringdist::stringdist(first(names(x)), last(names(x)), method = "hamming") )
  br_hm[is.infinite(br_hm)] <- 0
  br_hm_dist[[i]] <- br_hm
  
}

# Null.
br_hm_dist_null <- vector(mode = "list", length = length(stm_merg_null))
for (i in 1:length(stm_merg_null)) {
  
  br_hm_null <- sapply(stm_merg_null[[i]]$maps, function(x) stringdist::stringdist(first(names(x)), last(names(x)), method = "hamming") )
  br_hm_null[is.infinite(br_hm_null)] <- 0
  br_hm_dist_null[[i]] <- br_hm_null
  
}

# Get number of branches.
n_br = length(stm_amalg[[1]]$edge.length)

# Reorganize list of Hamming distances by branch.
# Data.
br_dist <- do.call(rbind, br_hm_dist)
colnames(br_dist) <- paste0("br_", 1:n_br)
br_dist <- as.list(as.data.frame(br_dist))

# Null.
br_dist_null <- do.call(rbind, br_hm_dist_null)
colnames(br_dist_null) <- paste0("br_", 1:n_br)
br_dist_null <- as.list(as.data.frame(br_dist_null))

# BOOTSTRAP SAMPLING #

# Set some parameters.
n_boot = 10000

# Data.
# Set a vector to store samples.
br_dist_boot <- numeric()

# Bootstrapping.
for (i in 1:n_boot) { br_dist_boot <- rbind(br_dist_boot, sapply(br_dist, function(x) mean(sample(x, length(x), replace = T)) )) }

# Save data.
save(br_dist, br_dist_boot, file = "boots.RDA")

# Null.
# Set a vector to store samples.
br_dist_boot_null <- numeric()

# Bootstrapping.
for (i in 1:n_boot) { br_dist_boot_null <- rbind(br_dist_boot_null, sapply(br_dist_null, function(x) mean(sample(x, length(x), replace = T)) )) }

# Save data.
save(br_dist_null, br_dist_boot_null, file = "boots_null.RDA")

# Calculate differences between Data and Null.
pval <- numeric()
h_diff <- vector(mode = "list", length = dim(br_dist_boot)[2])
for (i in 1:dim(br_dist_boot)[2]) {
  
  h_diff[[i]] <- br_dist_boot[,i] - br_dist_boot_null[,i]
  
  if (mean(h_diff[[i]]) < 0) {
    
    pval[[i]] <- mean(h_diff[[i]] >= 0)
    
  } else {
    
    pval[[i]] <- mean(h_diff[[i]] <= 0)
    
  }
    
}

# Calculate S statistics.
S_stat_adult <- sapply(h_diff, function(x) mean(x)/sqrt(var(x)) ) %>% round(.,2)

S_stat_adult[[5]] # Check S statistics for the focal edge.
pval[[5]] # Check p-val for the focal edge.
