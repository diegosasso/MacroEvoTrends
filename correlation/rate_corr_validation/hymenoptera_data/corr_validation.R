
# Load libraries.
library(ape)
library(phytools)
library(ontophylo)
library(dplyr)

# Import some functions.
source("R/helpers.R")

# Load Hymenoptera dataset.
load("data_hym/paramo_stm_final.RDA")

#--------------------------------------------------#

##############
# Check data #
##############

# Check tree.
pdf("hym_tree.pdf", height = 16, width = 12)
plot.phylo(hym_tree, show.tip.label = T, cex = 0.5, label.offset = 0.05)
edgelabels(frame = "none", col = "blue", cex = 0.8)
nodelabels(frame = "none", col = "red", cex = 0.8)
tiplabels(frame = "none", col = "purple", cex = 0.8)
dev.off()

# Check mean rates for individual traits from fitted models using empirical data.
rt_m <- sapply(fit_models_final, function(x) mean(rowSums(x, na.rm = T)) )
hist(rt_m, main = "", xlab = "Mean rates", 
     border = F, xlim = c(min(rt_m)*1.2, max(rt_m)*1.2), breaks = 50)
abline(v = mean(rt_m), lty = 2, col = "red")
abline(v = median(rt_m), lty = 2, col = "blue")
dev.off()

# Note: around 0.001

#--------------------------------------------------#

#########################
# Check empirical rates #
#########################

# Import data.
stm_amalg <- readRDS("data_hym/anat_ent.RDS")
head(stm_amalg)
length(stm_amalg$mesopectus)

# Merge tree bins.
stm_merg <- lapply(stm_amalg, function(x) merge_tree_cat_list(x) )

# Get max Hamming distances.
hm_dist <- lapply(stm_merg, function(x) hm_dist_br(x) )

# Calculate means across all stochastic maps.
hm_dist_m <- do.call(cbind, lapply(hm_dist, function(x) apply(do.call(cbind, x), 1, mean) ))

# Convert to rates.
hm_dist_rt <- hm_dist_m/hym_tree$edge.length

# Check summaries.
apply(hm_dist_rt, 2, min)
apply(hm_dist_rt, 2, max)
apply(hm_dist_rt, 2, median)

# Overall mean.
mean(apply(hm_dist_rt,1,mean))

# Note: around 0.03

#------- ST: new rate calculation
smaps <- stm_amalg$mesopectus

i=1
S <- smaps[[i]]
n.changes <- lapply(S$maps, function(x) length(x)-1) %>% unlist()


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

x <- get_branch_rate_across_smaps(stm_merg$mesopectus)
y <- hm_dist_rt[,'mesopectus']

plot(x,y)
#--------------------------------------------------#

# Set some parameters.
nstm = 100

#--------------------------------------------------#

####################
# Get branch rates #
####################

# Import data.
stm_amalg <- readRDS("data_hym/anat_ent.RDS")
stm_amalg$mesopectus

# Merge tree bins.
stm_merg <- lapply(stm_amalg, function(x) merge_tree_cat_list(x) )
stm_merg

# Get Hamming distances.
hm_dist <- lapply(stm_merg, function(x) hm_dist_br(x) )
hm_dist

# Calculate means across all stochastic maps.
hm_dist_m <- do.call(cbind, lapply(hm_dist, function(x) apply(do.call(cbind, x), 1, mean) ))
hm_dist_m
hm_dist_m

# Convert to rates.
hm_dist_rt <- lapply(hm_dist, function(x) lapply(x, function(y) y/hym_tree$edge.length ))
hm_dist_m_rt <- hm_dist_m/hym_tree$edge.length
hm_dist_m_rt

hm_dist_m_rt[,'mesopectus']
#--------------------------------------------------#

#####################
# PERMUTATION TESTS #
#####################

# Get branches leading to tips.
br_out <- which(!hym_tree$edge[,2] %in% c(1:length(hym_tree$tip.label)))
br_out

# Get number of partitions.
npart = length(colnames(hm_dist_m_rt))
npart

# Get all possible pairwise combinations.
perm_pairs <- combn(colnames(hm_dist_m_rt), 2)
perm_pairs

# Set matrix to store correlation values.
corr_vals <- matrix(NA, npart, npart)
colnames(corr_vals) <- rownames(corr_vals) <- colnames(hm_dist_m_rt)
corr_vals

# Set matrix to store p-values.
p_vals <- matrix(NA, npart, npart)
colnames(p_vals) <- rownames(p_vals) <- colnames(hm_dist_m_rt)
p_vals

# Initialize vectors to store values (to facilitate printing).
corr_vals_v <- setNames(vector(mode = "numeric", length = dim(perm_pairs)[2]), apply(perm_pairs, 2, function(z) paste0(z, collapse = " x ") ))
p_vals_v <- setNames(vector(mode = "numeric", length = dim(perm_pairs)[2]), apply(perm_pairs, 2, function(z) paste0(z, collapse = " x ") ))

# Set list to store permutation samples.
perm_samples <- vector(mode = "list", length = dim(perm_pairs)[2])
names(perm_samples) <- apply(perm_pairs, 2, function(z) paste0(z, collapse = " x ") )

i = 1
for (i in 1:dim(perm_pairs)[2]) {
  
  #####################
  # Extract variables #
  #####################

  # Raw changes.
  x <- hm_dist_m[, perm_pairs[1,i]]
  y <- hm_dist_m[, perm_pairs[2,i]]  
    
  # Rates.
  #x <- hm_dist_m_rt[, perm_pairs[1,i]]
  #y <- hm_dist_m_rt[, perm_pairs[2,i]]
  
  # Residuals. Raw changes.
  #x <- lm(hm_dist_m[, perm_pairs[1,i]] ~ hym_tree$edge.length)$residuals
  #y <- lm(hm_dist_m[, perm_pairs[2,i]] ~ hym_tree$edge.length)$residuals
  
  # Residuals. Rates.
  #x <- lm(hm_dist_m_rt[, perm_pairs[1,i]] ~ hym_tree$edge.length)$residuals
  #y <- lm(hm_dist_m_rt[, perm_pairs[2,i]] ~ hym_tree$edge.length)$residuals
  
  # Standardize variables.
  #x <- scale(x)
  #y <- scale(y)
  
  # Remove terminal branches.
  #x <- x[br_out]
  #y <- y[br_out]
  
  # Permutation test.
  pm_test <- perm_test(x, y, nperm = 10000, method = "spearman")
  
  # Store values in matrices.
  corr_vals[perm_pairs[1,i], perm_pairs[2,i]] <- pm_test$obs_corr
  p_vals[perm_pairs[1,i], perm_pairs[2,i]] <- pm_test$pval
  
  # Store values in vectors.
  corr_vals_v[[i]] <- pm_test$obs_corr
  p_vals_v[[i]] <- pm_test$pval
  
  # Store values in 
  perm_samples[[i]] <- pm_test$perm_corr
  
}

# Create directory to store data.
suppressWarnings(dir.create("data_out"))

# Change file name.
#f_name = "corr_raw"
#f_name = "corr_raw_no_tip"
f_name = "corr_rates"
#f_name = "corr_rates_no_tip"

#f_name = "corr_raw_res"
#f_name = "corr_raw_res_no_tip"
#f_name = "corr_rates_res"
#f_name = "corr_rates_res_no_tip"

# Save data.
save(hym_tree, stm_merg, hm_dist, hm_dist_m, hm_dist_rt, hm_dist_m_rt, 
     corr_vals, p_vals, corr_vals_v, p_vals_v, perm_pairs, file = paste0("data_out/", f_name, ".RDA"))

# Save CSV summary.
write.csv(cbind(round(corr_vals_v,2), p_vals_v), paste0(f_name, ".csv"))

#--------------------------------------------------#

# Notes:
#
# Good strategy:
#
# Summary (mean branch rates across all stochastic maps) + 
# standardizing (not really necessary since correlation is invariable to scale) +
# removing branches leading to tips since they can produce misleading results!
# Also, sometimes using residuals of the regression between raw changes and 
# branch lengths can help (but not always);

#--------------------------------------------------#
