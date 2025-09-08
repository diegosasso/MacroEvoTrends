
# Load libraries.
library(ape)
library(phytools)
library(ontophylo)

# Import some functions.
source("R/helpers.R")

# Load Hymenoptera dataset.
load("data_hym/paramo_stm_final.RDA")

#--------------------------------------------------#

# Import data.
stm_amalg <- readRDS("data_hym/anat_ent.RDS")

# Merge tree bins.
stm_merg <- lapply(stm_amalg, function(x) merge_tree_cat_list(x) )
stm_merg

# Get Hamming distances.
hm_dist <- lapply(stm_merg, function(x) hm_dist_br(x) )
hm_dist

# Calculate means across all stochastic maps.
hm_dist_m <- do.call(cbind, lapply(hm_dist, function(x) apply(do.call(cbind, x), 1, mean) ))
hm_dist_m

# Convert to rates.
hm_dist_rt <- lapply(hm_dist, function(x) lapply(x, function(y) y/hym_tree$edge.length ))
hm_dist_m_rt <- hm_dist_m/hym_tree$edge.length
hm_dist_m_rt

#--------------------------------------------------#

#########################
# Check possible biases #
#########################

# Plot branch rates vs. branch lengths (with and without outliers).

# Get an anatomical partition of interest.
# Raw changes (branch means across all stochastic maps).
br_anat <- hm_dist_m[,1]

# Branch rates (branch means across all stochastic maps).
br_anat <- hm_dist_m_rt[,1]

# Get non-outliers.
# Values below the top 5%.
#br_n_out <- which(br_anat <= quantile(br_anat, 0.95))

# Branches leading to tips.
br_n_out <- which(!hym_tree$edge[,2] %in% c(1:length(hym_tree$tip.label)))

# With outliers.
br_corr <- perm_test(br_anat, hym_tree$edge.length, nperm = 10000, method = "spearman")

# Check values.
br_corr$pval
br_corr$obs_corr

# Scatter plot.
plot(br_anat, hym_tree$edge.length, xlab = "Br_changes/Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_corr$perm_corr, main = "Br_changes/Rates x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round((abs(br_corr$obs_corr)*1.5),2),round((abs(br_corr$obs_corr)*1.5),3)), breaks = 50)
abline(v = br_corr$obs_corr, lty = 2, col = "red")

#----------#

# Without outliers.
br_corr_c <- perm_test(br_anat[br_n_out], hym_tree$edge.length[br_n_out], nperm = 10000, method = "spearman")

# Check values.
br_corr_c$pval
br_corr_c$obs_corr

# Scatter plot.
plot(br_anat[br_n_out], hym_tree$edge.length[br_n_out], xlab = "Br_changes/Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_corr_c$perm_corr, main = "Br_changes/Rates x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round((abs(br_corr_c$obs_corr)*1.5),2),round((abs(br_corr_c$obs_corr)*1.5),3)), breaks = 50)
abline(v = br_corr_c$obs_corr, lty = 2, col = "red")

#--------------------------------------------------#

# NOTES:
#
# By removing the top 5% outliers (likely the same/close branches purposely simulated with high rates)
# the residual branch rates clearly show correlation with branch lengths, as expected.
# This applies when assessing "good" traits + "noisy" traits and only "noisy" traits;
# Thus it seems to be an effect of branch lengths on rates;
# Removing branches leading to tips seems to improve results since such branches
# are usually reconstructed with very low rates from empirical data given that
# autapomorphies are usually not scored;

#--------------------------------------------------#
