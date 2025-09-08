
# Load libraries.
library(ape)
library(phytools)
library(ontophylo)

# Import some functions.
source("R/utils.R")
source("R/helpers.R")

# Import tree.
tree <- readRDS("tree100.RDS")

# Rename tip labels to make it easier to visualize.
tree$tip.label <- paste0("t", 1:length(tree$tip.label))

#--------------------------------------------------#
# CONDITIONED ON DATA
#--------------------------------------------------#

# Import data.
stm_amalg <- readRDS("data_out/stm_amalg.RDS")

# Merge tree bins.
stm_merg <- lapply(stm_amalg, function(x) merge_tree_cat_list(x) )

# Get max Hamming distances.
hm_dist <- lapply(stm_merg, function(x) hm_dist_br(x) )

# Calculate means across all stochastic maps.
hm_dist_m <- do.call(cbind, lapply(hm_dist, function(x) apply(do.call(cbind, x), 1, mean) ))

# Convert to rates.
hm_dist <- lapply(hm_dist, function(x) lapply(x, function(y) y/tree$edge.length ))
hm_dist_m <- hm_dist_m/tree$edge.length

#--------------------------------------------------#
# UNCONDITIONED ON DATA
#--------------------------------------------------#

# Import data.
stm_un_amalg <- readRDS("data_out/stm_amalg_uncond.RDS")

# Merge tree bins.
stm_un_merg <- lapply(stm_un_amalg, function(x) merge_tree_cat_list(x) )

# Get max Hamming distances.
hm_dist_un <- lapply(stm_un_merg, function(x) hm_dist_br(x) )

# Calculate means across all stochastic maps.
hm_dist_un_m <- do.call(cbind, lapply(hm_dist_un, function(x) apply(do.call(cbind, x), 1, mean) ))

# Convert to rates.
hm_dist_un <- lapply(hm_dist_un, function(x) lapply(x, function(y) y/tree$edge.length ))
hm_dist_un_m <- hm_dist_un_m/tree$edge.length

#--------------------------------------------------#
# INVESTIGATE POSSIBLE BIASES
#--------------------------------------------------#

## BRANCH LENGTHS #

# CONDITIONED #
# Scatter plots.
plot(hm_dist_m[,"A1"], tree$edge.length, xlab = "A1", ylab = "Brlen", pch = 19, col = "red")
dev.off()

plot(hm_dist_m[,"A2"], tree$edge.length, xlab = "A1", ylab = "Brlen", pch = 19, col = "red")
dev.off()

plot(hm_dist_m[,"B1"], tree$edge.length, xlab = "A1", ylab = "Brlen", pch = 19, col = "red")
dev.off()

plot(hm_dist_m[,"B2"], tree$edge.length, xlab = "A1", ylab = "Brlen", pch = 19, col = "red")
dev.off()

# UNCONDITIONED #
# Scatter plots.
plot(hm_dist_un_m[,"A1"], tree$edge.length, xlab = "A1", ylab = "Brlen", pch = 19, col = "red")
dev.off()

plot(hm_dist_un_m[,"A2"], tree$edge.length, xlab = "A1", ylab = "Brlen", pch = 19, col = "red")
dev.off()

plot(hm_dist_un_m[,"B1"], tree$edge.length, xlab = "A1", ylab = "Brlen", pch = 19, col = "red")
dev.off()

plot(hm_dist_un_m[,"B2"], tree$edge.length, xlab = "A1", ylab = "Brlen", pch = 19, col = "red")
dev.off()

# Notes: yes, obviously there are more changes onto longer branches!

#--------------------------------------------------#
## TRY CORRECTIONS ##
#--------------------------------------------------#

# Get only difference in changes
#hm_dist_c_m <- sqrt((hm_dist_m - hm_dist_un_m)^2)
#hm_dist_c_m <- hm_dist_m - hm_dist_un_m
hm_dist_c_m <- hm_dist_m/hm_dist_un_m
hm_dist_c_m[is.nan(hm_dist_c_m)] <- 0

#--------------------------------------------------#
# CHECK FOR CORRELATIONS #
#--------------------------------------------------#

############################
### PERMUTATION APPROACH ###
############################

#--------------------------------------------------#

# SHOULD BE CORRELATED

#--------------------------------------------------#

##################
# Groups A1 x A2 #
##################

# Permutation test.
perm_t_GA1x2_c <- perm_test(hm_dist_c_m[,"A1"], hm_dist_c_m[,"A2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1x2_c$pval

# Observed correlation.
perm_t_GA1x2_c$obs_corr

# Scatter plot.
plot(hm_dist_c_m[,"A1"], hm_dist_c_m[,"A2"], xlab = "A1", ylab = "A2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1x2_c$perm_corr, main = "A1 x A2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1x2_c$obs_corr)*1.3),2),round(((perm_t_GA1x2_c$obs_corr)*1.3),2)), breaks = 50)
abline(v = perm_t_GA1x2_c$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups B1 x B2 #
##################

# Permutation test.
perm_t_GB1x2_c <- perm_test(hm_dist_c_m[,"B1"], hm_dist_c_m[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GB1x2_c$pval

# Observed correlation.
perm_t_GB1x2_c$obs_corr

# Scatter plot.
plot(hm_dist_c_m[,"B1"], hm_dist_c_m[,"B2"], xlab = "B1", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GB1x2_c$perm_corr, main = "B1 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GB1x2_c$obs_corr)*1.3),2),round(((perm_t_GB1x2_c$obs_corr)*1.3),2)), breaks = 50)
abline(v = perm_t_GB1x2_c$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------------------------------------#

# SHOULD BE NOT CORRELATED!!!

#--------------------------------------------------#

##################
# Groups A1 x B1 #
##################

# Permutation test.
perm_t_GA1xB1_c <- perm_test(hm_dist_c_m[,"A1"], hm_dist_c_m[,"B1"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1xB1_c$pval

# Observed correlation.
perm_t_GA1xB1_c$obs_corr

# Scatter plot.
plot(hm_dist_c_m[,"A1"], hm_dist_c_m[,"B1"], xlab = "A1", ylab = "B1", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1xB1_c$perm_corr, main = "A1 x B1", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1xB1_c$obs_corr)*1.3),2),round(((perm_t_GA1xB1_c$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA1xB1_c$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups A1 x B2 #
##################

# Permutation test.
perm_t_GA1xB2_c <- perm_test(hm_dist_c_m[,"A1"], hm_dist_c_m[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1xB2_c$pval

# Observed correlation.
perm_t_GA1xB2_c$obs_corr

# Scatter plot.
plot(hm_dist_c_m[,"A1"], hm_dist_c_m[,"B2"], xlab = "A1", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1xB2_c$perm_corr, main = "A1 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1xB2_c$obs_corr)*1.3),2),round(((perm_t_GA1xB2_c$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA1xB2_c$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups A2 x B1 #
##################

# Permutation test.
perm_t_GA2xB1_c <- perm_test(hm_dist_c_m[,"A2"], hm_dist_c_m[,"B1"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA2xB1_c$pval

# Observed correlation.
perm_t_GA2xB1_c$obs_corr

# Scatter plot.
plot(hm_dist_c_m[,"A2"], hm_dist_c_m[,"B1"], xlab = "A2", ylab = "B1", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA2xB1_c$perm_corr, main = "A2 x B1", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA2xB1_c$obs_corr)*1.3),2),round(((perm_t_GA2xB1_c$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA2xB1_c$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups A2 x B2 #
##################

# Permutation test.
perm_t_GA2xB2_c <- perm_test(hm_dist_c_m[,"A2"], hm_dist_c_m[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA2xB2_c$pval

# Observed correlation.
perm_t_GA2xB2_c$obs_corr

# Scatter plot.
plot(hm_dist_c_m[,"A2"], hm_dist_c_m[,"B2"], xlab = "A2", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA2xB2_c$perm_corr, main = "A2 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA2xB2_c$obs_corr)*1.3),2),round(((perm_t_GA2xB2_c$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA2xB2_c$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------------------------------------#

# Save data.
save(perm_t_GA1x2_c, perm_t_GB1x2_c, perm_t_GA1xB1_c, perm_t_GA1xB2_c, perm_t_GA2xB1_c, perm_t_GA2xB2_c,
     file = "data_out/permut_tests_corrected.RDA")

#--------------------------------------------------#

#########################
### MORLON'S APPROACH ###
#########################

# Set list to re-scaled trees.
tree_set_c <- setNames(vector(mode = "list", length = dim(hm_dist_c_m)[2]), colnames(hm_dist_c_m))
tree_set_c

j = 1
# Loop over.
for (j in 1:length(tree_set_c)) {
  
  # Copy original tree.
  tree_c <- tree
  
  # Re-scale branches by mean observed changes. Replace branch lengths from original tree.
  tree_c$edge.length <- hm_dist_c_m[,j]
  #tree_c$edge.length[tree_c$edge.length == 0] <- 1e-3
  
  # Store re-scaled tree.
  tree_set_c[[j]] <- tree_c
  
}

# Load packages.
library(RPANDA)

# Calculate JSD among trees. Probably a bad idea with too few trees!
jsd_c = JSDtree(tree_set_c, meth = c("standard"))

# Get clusters.
JSDtree_cluster(jsd_c)
