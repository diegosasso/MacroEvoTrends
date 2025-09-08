
# Load libraries.
library(ape)
library(phytools)
library(rphenoscate)
library(ontophylo)
library(corHMM)

# Import some functions.
source("R/utils.R")
source("R/helpers.R")

# Load Hymenoptera dataset.
load("data_hym/paramo_stm_final.RDA")

# Import data.
stm_amalg <- readRDS("data_out/stm_amalg.RDS")
stm_amalg <- readRDS("data_out/stm_amalg_alt1.RDS")
stm_amalg <- readRDS("data_out/stm_amalg_alt2.RDS")

#--------------------------------------------------#

####################
# Get branch rates #
####################

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

# Use raw changes instead (just reusing object to not change code below).
#hm_dist_m_rt <- hm_dist_m

#--------------------------------------------------#

###################################
# Check possible bias/corrections #
###################################

# Plot branch rates vs. branch lengths (with and without top 5% outliers).

######
# A1 #
######

# Get non-outliers.
#br_out1 <- which(hm_dist_m_rt[,"A1"] <= quantile(hm_dist_m_rt[,"A1"], 0.95))
br_out1 <- which(!hym_tree$edge[,2] %in% c(1:length(hym_tree$tip.label)))

# With top 5% outliers.
br_cor_a1 <- perm_test(hm_dist_m_rt[,"A1"], hym_tree$edge.length, nperm = 10000, method = "spearman")

# Check values.
br_cor_a1$pval
br_cor_a1$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"A1"], hym_tree$edge.length, xlab = "Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_cor_a1$perm_corr, main = "A1 x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round(((br_cor_a1$obs_corr)*5.0),2),round(((br_cor_a1$obs_corr)*5.0),3)), breaks = 50)
abline(v = br_cor_a1$obs_corr, lty = 2, col = "red")

#----------#

# Without top 5% outliers.
br_cor_a1_c <- perm_test(hm_dist_m_rt[,"A1"][br_out1], hym_tree$edge.length[br_out1], nperm = 10000, method = "spearman")

# Check values.
br_cor_a1_c$pval
br_cor_a1_c$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"A1"][br_out1], hym_tree$edge.length[br_out1], xlab = "Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_cor_a1_c$perm_corr, main = "A1 x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round(((br_cor_a1_c$obs_corr)*3.0),2),round(((br_cor_a1_c$obs_corr)*3.0),3)), breaks = 50)
abline(v = br_cor_a1_c$obs_corr, lty = 2, col = "red")

#--------------------#

######
# A2 #
######

# Get non-outliers.
#br_out2 <- which(hm_dist_m_rt[,"A2"] <= quantile(hm_dist_m_rt[,"A2"], 0.95))
br_out2 <- which(!hym_tree$edge[,2] %in% c(1:length(hym_tree$tip.label)))

# With top 5% outliers.
br_cor_a2 <- perm_test(hm_dist_m_rt[,"A2"], hym_tree$edge.length, nperm = 10000, method = "spearman")

# Check values.
br_cor_a2$pval
br_cor_a2$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"A2"], hym_tree$edge.length, xlab = "Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_cor_a2$perm_corr, main = "A2 x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round(((br_cor_a2$obs_corr)*5.0),2),round(((br_cor_a2$obs_corr)*5.0),3)), breaks = 50)
abline(v = br_cor_a2$obs_corr, lty = 2, col = "red")

#----------#

# Without top 5% outliers.
br_cor_a2_c <- perm_test(hm_dist_m_rt[,"A2"][br_out2], hym_tree$edge.length[br_out2], nperm = 10000, method = "spearman")

# Check values.
br_cor_a2_c$pval
br_cor_a2_c$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"A2"][br_out2], hym_tree$edge.length[br_out2], xlab = "Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_cor_a2_c$perm_corr, main = "A2 x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round(((br_cor_a2_c$obs_corr)*3.0),2),round(((br_cor_a2_c$obs_corr)*3.0),3)), breaks = 50)
abline(v = br_cor_a2_c$obs_corr, lty = 2, col = "red")

#--------------------#

######
# B1 #
######

# Get non-outliers.
#br_out3 <- which(hm_dist_m_rt[,"B1"] <= quantile(hm_dist_m_rt[,"B1"], 0.95))
br_out3 <- which(!hym_tree$edge[,2] %in% c(1:length(hym_tree$tip.label)))

# With top 5% outliers.
br_cor_b1 <- perm_test(hm_dist_m_rt[,"B1"], hym_tree$edge.length, nperm = 10000, method = "spearman")

# Check values.
br_cor_b1$pval
br_cor_b1$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"B1"], hym_tree$edge.length, xlab = "Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_cor_b1$perm_corr, main = "B1 x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round(((br_cor_b1$obs_corr)*5.0),2),round(((br_cor_b1$obs_corr)*5.0),3)), breaks = 50)
abline(v = br_cor_b1$obs_corr, lty = 2, col = "red")

#----------#

# Without top 5% outliers.
br_cor_b1_c <- perm_test(hm_dist_m_rt[,"B1"][br_out3], hym_tree$edge.length[br_out3], nperm = 10000, method = "spearman")

# Check values.
br_cor_b1_c$pval
br_cor_b1_c$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"B1"][br_out3], hym_tree$edge.length[br_out3], xlab = "Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_cor_b1_c$perm_corr, main = "B1 x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round(((br_cor_b1_c$obs_corr)*3.0),2),round(((br_cor_b1_c$obs_corr)*3.0),3)), breaks = 50)
abline(v = br_cor_b1_c$obs_corr, lty = 2, col = "red")

#--------------------#

######
# B2 #
######

# Get non-outliers.
#br_out4 <- which(hm_dist_m_rt[,"B2"] <= quantile(hm_dist_m_rt[,"B2"], 0.95))
br_out4 <- which(!hym_tree$edge[,2] %in% c(1:length(hym_tree$tip.label)))

# With top 5% outliers.
br_cor_b2 <- perm_test(hm_dist_m_rt[,"B2"], hym_tree$edge.length, nperm = 10000, method = "spearman")

# Check values.
br_cor_b2$pval
br_cor_b2$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"B2"], hym_tree$edge.length, xlab = "Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_cor_b2$perm_corr, main = "B2 x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round(((br_cor_b2$obs_corr)*5.0),2),round(((br_cor_b2$obs_corr)*5.0),3)), breaks = 50)
abline(v = br_cor_b2$obs_corr, lty = 2, col = "red")

#----------#

# Without top 5% outliers.
br_cor_b2_c <- perm_test(hm_dist_m_rt[,"B2"][br_out4], hym_tree$edge.length[br_out4], nperm = 10000, method = "spearman")

# Check values.
br_cor_b2_c$pval
br_cor_b2_c$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"B2"][br_out4], hym_tree$edge.length[br_out4], xlab = "Rates", ylab = "Brlen", pch = 19, col = "purple")

# Histogram.
hist(br_cor_b2_c$perm_corr, main = "B2 x Brlen", xlab = "correlation", 
     border = F, xlim = c(-round(((br_cor_b2_c$obs_corr)*3.0),2),round(((br_cor_b2_c$obs_corr)*3.0),3)), breaks = 50)
abline(v = br_cor_b2_c$obs_corr, lty = 2, col = "red")

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



