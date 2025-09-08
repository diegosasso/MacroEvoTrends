
###############################################################
### SAME AS QUICK TEST 8 BUT USING RATES INSTEAD OF CHANGES ###
###############################################################

# Load libraries.
library(ape)
library(phytools)
library(rphenoscate)
library(ontophylo)
library(corHMM)

# Import some functions.
source("correlation/rate_corr_validation/R/utils.R")
source("correlation/rate_corr_validation/R/helpers.R")

# Create folders to store data.
suppressWarnings(dir.create("data_out"))

#--------------------------------------------------#

# Import tree.
tree <- readRDS("correlation/rate_corr_validation/tree100.RDS")
plot(tree)

# Rename tip labels to make it easier to visualize.
tree$tip.label <- paste0("t", 1:length(tree$tip.label))

# Check tree.
#pdf("tree.pdf", height = 16, width = 12)
plot.phylo(tree, show.tip.label = T, cex = 0.5, label.offset = 0.05)
edgelabels(frame = "none", col = "blue", cex = 0.8)
nodelabels(frame = "none", col = "red", cex = 0.8)
tiplabels(frame = "none", col = "purple", cex = 0.8)
dev.off()

#--------------------------------------------------#

#------------------------------#
# Set model.
#------------------------------#

# Set base Q-matrix.
Q_base <- initQ(c(0, 1), c(0.01,0.01))
Q_base

# Set Q-matrix of shift.
# Note that Q_jump might be reordered depending on initial state.
Q_jump <- initQ(c(0, 1), c(1, 1e-10))
Q_jump

#--------------------------------------------------#

# Set some parameters.
nsim = 10
nstm = 100

#--------------------------------------------------#

#------------------------------#
# Simulate data.
#------------------------------#

# Define a set of focal branches for shifts.
br = c(1, 2, 3, 90, 92, 150)

# Set a list to store simulated tip data.
sim_data <- setNames(vector(mode = "list", length = length(br)), paste0("BR_", br))

# Set a list to store sampled character histories.
stm_data <- setNames(vector(mode = "list", length = length(br)), paste0("BR_", br))

# Simulate data and sample character histories.
for (i in 1:length(br)) {
  
  # Simulate tip data.
  sim_data[[i]] <- sim_info_data(tree = tree, focal_edge = br[[i]], Q_base = Q_base, Q_jump = Q_jump, nsim = nsim)
  
  # Sample character histories conditioned on tip data.
  stm_data[[i]] <- sample_cond(tree, sim_data[[i]], nsim, nstm)
  
}

# Save data.
save(sim_data, stm_data, file = "data_out/test_data.RDA")

# Check some maps.
# BR_1.
densityMap(stm_data$BR_1$smaps$C1)
densityMap(stm_data$BR_1$smaps$C2)
densityMap(stm_data$BR_1$smaps$C3)

# BR_2.
densityMap(stm_data$BR_2$smaps$C1)
densityMap(stm_data$BR_2$smaps$C2)
densityMap(stm_data$BR_2$smaps$C3)

# BR_3.
densityMap(stm_data$BR_3$smaps$C1)
densityMap(stm_data$BR_3$smaps$C2)
densityMap(stm_data$BR_3$smaps$C3)

# BR_90.
densityMap(stm_data$BR_90$smaps$C1)
densityMap(stm_data$BR_90$smaps$C2)
densityMap(stm_data$BR_90$smaps$C3)

# BR_92.
densityMap(stm_data$BR_92$smaps$C1)
densityMap(stm_data$BR_92$smaps$C2)
densityMap(stm_data$BR_92$smaps$C3)

# BR_150.
densityMap(stm_data$BR_150$smaps$C1)
densityMap(stm_data$BR_150$smaps$C2)
densityMap(stm_data$BR_150$smaps$C3)

#--------------------------------------------------#

#------------------------------#
# Split datasets.
#------------------------------#

# Split datasets in half to emulate two groups of traits highly correlated.
# EXPLANATION NOTES:
# the first group (A) comprise binary traits changing along branches 1, 2, and 3;
# 10 binary traits simulated for each case (i.e., focal branch), thus 30 in total.
# this emulates a sequence of branches with high evolutionary rates leading to 
# the clade comprising tips t1 to t 15.

# the second group (B) comprise binary traits changing along branches 90, 92, and 150;
# 10 binary traits simulated for each case (i.e., focal branch), thus 30 in total.
# this emulates a sequence of branches with high evolutionary rates leading to 
# the clade comprising tips t46 to t 100.

# Then each group (A) and (B) were split in half, A1-A2 and B1-B2 (15 traits each).
# PARAMO pipeline is applied and characters amalgamated resulting in 4 trees (100 stochastic maps each).
# Branch changes were extracted as a proxy for branch rates.
# Expectation:
# (1) groups A1 and A2 should be highly correlated; same as groups B1 and B2;
# (2) groups A1 x B1 or B2 and A2 x B1 or B2 should not be correlated (i.e., no branches with high rates in common);

# Group A1.
stm_GA1 <- lapply(stm_data[c(1,2,3)], function(x) x$smaps[c(1:5)] )

# Group A2.
stm_GA2 <- lapply(stm_data[c(1,2,3)], function(x) x$smaps[c(6:10)] )

# Group B1.
stm_GB1 <- lapply(stm_data[c(4,5,6)], function(x) x$smaps[c(1:5)] )

# Group B2.
stm_GB2 <- lapply(stm_data[c(4,5,6)], function(x) x$smaps[c(6:10)] )

#--------------------------------------------------#

#------------------------------#
# Amalgamate characters.
#------------------------------#

# Set parameters.
res = 200

# Consolidate datasets and rename characters.
stm_GA1 <- setNames(do.call(c, stm_GA1), paste0("C", 1:sum(sapply(stm_GA1, length))))
stm_GA2 <- setNames(do.call(c, stm_GA2), paste0("C", 1:sum(sapply(stm_GA2, length))))
stm_GB1 <- setNames(do.call(c, stm_GB1), paste0("C", 1:sum(sapply(stm_GB1, length))))
stm_GB2 <- setNames(do.call(c, stm_GB2), paste0("C", 1:sum(sapply(stm_GB2, length))))

# Discretize maps.
stm_GA1_discr <- lapply(stm_GA1, function(x) discr_Simmap_all(x, res) )
stm_GA2_discr <- lapply(stm_GA2, function(x) discr_Simmap_all(x, res) )
stm_GB1_discr <- lapply(stm_GB1, function(x) discr_Simmap_all(x, res) )
stm_GB2_discr <- lapply(stm_GB2, function(x) discr_Simmap_all(x, res) )

# Amalgamate maps.
stm_GA1_amalg <- paramo.list(names(stm_GA1_discr), tree.list = stm_GA1_discr, ntrees = nstm)
stm_GA2_amalg <- paramo.list(names(stm_GA2_discr), tree.list = stm_GA2_discr, ntrees = nstm)
stm_GB1_amalg <- paramo.list(names(stm_GB1_discr), tree.list = stm_GB1_discr, ntrees = nstm)
stm_GB2_amalg <- paramo.list(names(stm_GB2_discr), tree.list = stm_GB2_discr, ntrees = nstm)

# Consolidate datasets.
stm_amalg <- list("A1" = stm_GA1_amalg,
                  "A2" = stm_GA2_amalg,
                  "B1" = stm_GB1_amalg,
                  "B2" = stm_GB2_amalg)

# Save data.
#saveRDS(stm_amalg, "correlation/rate_corr_validation/data_out/stm_amalg.RDS")

#--------------------------------------------------#

#------------------------------#
# Check for correlations.
#------------------------------#

# Import data.
stm_amalg <- readRDS("data_out/stm_amalg.RDS")

# Merge tree bins.
stm_merg <- lapply(stm_amalg, function(x) merge_tree_cat_list(x) )
stm_merg

# Get max Hamming distances.
hm_dist <- lapply(stm_merg, function(x) hm_dist_br(x) )
hm_dist

# Calculate means across all stochastic maps.
hm_dist_m <- do.call(cbind, lapply(hm_dist, function(x) apply(do.call(cbind, x), 1, mean) ))
hm_dist_m

# Convert to rates.
hm_dist <- lapply(hm_dist, function(x) lapply(x, function(y) y/tree$edge.length ))
hm_dist_m <- hm_dist_m/tree$edge.length

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
perm_t_GA1x2 <- perm_test(hm_dist_m[,"A1"], hm_dist_m[,"A2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1x2$pval

# Observed correlation.
perm_t_GA1x2$obs_corr

# Scatter plot.
plot(hm_dist_m[,"A1"], hm_dist_m[,"A2"], xlab = "A1", ylab = "A2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1x2$perm_corr, main = "A1 x A2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1x2$obs_corr)*1.3),2),round(((perm_t_GA1x2$obs_corr)*1.3),2)), breaks = 50)
abline(v = perm_t_GA1x2$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups B1 x B2 #
##################

# Permutation test.
perm_t_GB1x2 <- perm_test(hm_dist_m[,"B1"], hm_dist_m[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GB1x2$pval

# Observed correlation.
perm_t_GB1x2$obs_corr

# Scatter plot.
plot(hm_dist_m[,"B1"], hm_dist_m[,"B2"], xlab = "B1", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GB1x2$perm_corr, main = "B1 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GB1x2$obs_corr)*1.3),2),round(((perm_t_GB1x2$obs_corr)*1.3),2)), breaks = 50)
abline(v = perm_t_GB1x2$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------------------------------------#

# SHOULD BE NOT CORRELATED!!!

#--------------------------------------------------#

##################
# Groups A1 x B1 #
##################

# Permutation test.
perm_t_GA1xB1 <- perm_test(hm_dist_m[,"A1"], hm_dist_m[,"B1"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1xB1$pval

# Observed correlation.
perm_t_GA1xB1$obs_corr

# Scatter plot.
plot(hm_dist_m[,"A1"], hm_dist_m[,"B1"], xlab = "A1", ylab = "B1", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1xB1$perm_corr, main = "A1 x B1", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1xB1$obs_corr)*1.3),2),round(((perm_t_GA1xB1$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA1xB1$obs_corr, lty = 2, col = "red")
dev.off()

#----------#
##################
# Groups A1 x B1 # perm_test_Double
##################

# Permutation test.
perm_t_GA1xB1_Double <- perm_test_Double(hm_dist_m[,"A1"], hm_dist_m[,"B1"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1xB1_Double$pval

##################
# Groups A1 x B2 #
##################

# Permutation test.
perm_t_GA1xB2 <- perm_test(hm_dist_m[,"A1"], hm_dist_m[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1xB2$pval

# Observed correlation.
perm_t_GA1xB2$obs_corr

# Scatter plot.
plot(hm_dist_m[,"A1"], hm_dist_m[,"B2"], xlab = "A1", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1xB2$perm_corr, main = "A1 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1xB2$obs_corr)*1.3),2),round(((perm_t_GA1xB2$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA1xB2$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups A2 x B1 #
##################

# Permutation test.
perm_t_GA2xB1 <- perm_test(hm_dist_m[,"A2"], hm_dist_m[,"B1"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA2xB1$pval

# Observed correlation.
perm_t_GA2xB1$obs_corr

# Scatter plot.
plot(hm_dist_m[,"A2"], hm_dist_m[,"B1"], xlab = "A2", ylab = "B1", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA2xB1$perm_corr, main = "A2 x B1", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA2xB1$obs_corr)*1.3),2),round(((perm_t_GA2xB1$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA2xB1$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups A2 x B2 #
##################

# Permutation test.
perm_t_GA2xB2 <- perm_test(hm_dist_m[,"A2"], hm_dist_m[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA2xB2$pval

# Observed correlation.
perm_t_GA2xB2$obs_corr

# Scatter plot.
plot(hm_dist_m[,"A2"], hm_dist_m[,"B2"], xlab = "A2", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA2xB2$perm_corr, main = "A2 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA2xB2$obs_corr)*1.3),2),round(((perm_t_GA2xB2$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA2xB2$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------------------------------------#

# Save data.
save(perm_t_GA1x2, perm_t_GB1x2, perm_t_GA1xB1, perm_t_GA1xB2, perm_t_GA2xB1, perm_t_GA2xB2,
     file = "data_out/permut_tests.RDA")

#--------------------------------------------------#

#########################
### MORLON'S APPROACH ###
#########################

# Set list to re-scaled trees.
tree_set <- setNames(vector(mode = "list", length = dim(hm_dist_m)[2]), colnames(hm_dist_m))
tree_set

j = 1
# Loop over.
for (j in 1:length(tree_set)) {
  
  # Copy original tree.
  tree_c <- tree
  
  # Re-scale branches by mean observed changes. Replace branch lengths from original tree.
  tree_c$edge.length <- hm_dist_m[,j]
  
  # Store re-scaled tree.
  tree_set[[j]] <- tree_c
  tree_c$edge.length[tree_c$edge.length == 0] <- 1e-3
  
}

# Load packages.
library(RPANDA)

# Calculate JSD among trees. Probably a bad idea with too few trees!
jsd = JSDtree(tree_set, meth = c("standard"))

# Get clusters.
JSDtree_cluster(jsd)

# Notes:
# Probably not working properly since several branches have 0 length.



################################################################################



#--------------------------------------------------#

# CHECK CORRELATIONS WITH UNCONDITIONED SAMPLES #

#--------------------------------------------------#

# Load data.
load("data_out/test_data.RDA")


# Set some parameters.
nsim = 10
nstm = 100
res = 200

#------------------------------#

####################
# Re-organize data #
####################

# Get fitted models.
# Group A1.
models_GA1 <- lapply(stm_data[c(1,2,3)], function(x) x$models[c(1:5)] )

# Group A2.
models_GA2 <- lapply(stm_data[c(1,2,3)], function(x) x$models[c(6:10)] )

# Group B1.
models_GB1 <- lapply(stm_data[c(4,5,6)], function(x) x$models[c(1:5)] )

# Group B2.
models_GB2 <- lapply(stm_data[c(4,5,6)], function(x) x$models[c(6:10)] )

#------------------------------#

# Sample character histories unconditioned on data.
stm_GA1_un <- lapply(models_GA1, function(x) sim_uncond(tree, fitted_models = x, nsim = (nsim/2), nstm = nstm))
stm_GA2_un <- lapply(models_GA2, function(x) sim_uncond(tree, fitted_models = x, nsim = (nsim/2), nstm = nstm))
stm_GB1_un <- lapply(models_GB1, function(x) sim_uncond(tree, fitted_models = x, nsim = (nsim/2), nstm = nstm))
stm_GB2_un <- lapply(models_GB2, function(x) sim_uncond(tree, fitted_models = x, nsim = (nsim/2), nstm = nstm))

# Consolidate datasets and rename characters.
stm_GA1_un <- setNames(do.call(c, stm_GA1_un), paste0("C", 1:sum(sapply(stm_GA1_un, length))))
stm_GA2_un <- setNames(do.call(c, stm_GA2_un), paste0("C", 1:sum(sapply(stm_GA2_un, length))))
stm_GB1_un <- setNames(do.call(c, stm_GB1_un), paste0("C", 1:sum(sapply(stm_GB1_un, length))))
stm_GB2_un <- setNames(do.call(c, stm_GB2_un), paste0("C", 1:sum(sapply(stm_GB2_un, length))))

# Save data.
save(stm_GA1_un, stm_GA2_un, stm_GB1_un, stm_GB2_un, file = "data_out/test_data_uncond.RDA")

#------------------------------#

# Discretize maps.
stm_GA1_un_discr <- lapply(stm_GA1_un, function(x) discr_Simmap_all(x, res) )
stm_GA2_un_discr <- lapply(stm_GA2_un, function(x) discr_Simmap_all(x, res) )
stm_GB1_un_discr <- lapply(stm_GB1_un, function(x) discr_Simmap_all(x, res) )
stm_GB2_un_discr <- lapply(stm_GB2_un, function(x) discr_Simmap_all(x, res) )

# Amalgamate maps.
stm_GA1_un_amalg <- paramo.list(names(stm_GA1_un_discr), tree.list = stm_GA1_un_discr, ntrees = nstm)
stm_GA2_un_amalg <- paramo.list(names(stm_GA2_un_discr), tree.list = stm_GA2_un_discr, ntrees = nstm)
stm_GB1_un_amalg <- paramo.list(names(stm_GB1_un_discr), tree.list = stm_GB1_un_discr, ntrees = nstm)
stm_GB2_un_amalg <- paramo.list(names(stm_GB2_un_discr), tree.list = stm_GB2_un_discr, ntrees = nstm)

# Consolidate datasets.
stm_un_amalg <- list("A1" = stm_GA1_un_amalg,
                     "A2" = stm_GA2_un_amalg,
                     "B1" = stm_GB1_un_amalg,
                     "B2" = stm_GB2_un_amalg)

# Save data.
saveRDS(stm_un_amalg, "data_out/stm_amalg_uncond.RDS")

#--------------------------------------------------#

#------------------------------#
# Check for correlations.
#------------------------------#

# Import data.
stm_un_amalg <- readRDS("data_out/stm_amalg_uncond.RDS")

# Merge tree bins.
stm_un_merg <- lapply(stm_un_amalg, function(x) merge_tree_cat_list(x) )
stm_un_merg

# Get max Hamming distances.
hm_dist_un <- lapply(stm_un_merg, function(x) hm_dist_br(x) )
hm_dist_un

# Calculate means across all stochastic maps.
hm_dist_un_m <- do.call(cbind, lapply(hm_dist_un, function(x) apply(do.call(cbind, x), 1, mean) ))
hm_dist_un_m

# Convert to rates.
hm_dist_un <- lapply(hm_dist_un, function(x) lapply(x, function(y) y/tree$edge.length ))
hm_dist_un_m <- hm_dist_un_m/tree$edge.length

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
perm_t_GA1x2_un <- perm_test(hm_dist_un_m[,"A1"], hm_dist_un_m[,"A2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1x2_un$pval

# Observed correlation.
perm_t_GA1x2_un$obs_corr

# Scatter plot.
plot(hm_dist_un_m[,"A1"], hm_dist_un_m[,"A2"], xlab = "A1", ylab = "A2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1x2_un$perm_corr, main = "A1 x A2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1x2_un$obs_corr)*1.3),2),round(((perm_t_GA1x2_un$obs_corr)*1.3),2)), breaks = 50)
abline(v = perm_t_GA1x2_un$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups B1 x B2 #
##################

# Permutation test.
perm_t_GB1x2_un <- perm_test(hm_dist_un_m[,"B1"], hm_dist_un_m[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GB1x2_un$pval

# Observed correlation.
perm_t_GB1x2_un$obs_corr

# Scatter plot.
plot(hm_dist_un_m[,"B1"], hm_dist_un_m[,"B2"], xlab = "B1", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GB1x2_un$perm_corr, main = "B1 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GB1x2_un$obs_corr)*1.3),2),round(((perm_t_GB1x2_un$obs_corr)*1.3),2)), breaks = 50)
abline(v = perm_t_GB1x2_un$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------------------------------------#

# SHOULD BE NOT CORRELATED!!!

#--------------------------------------------------#

##################
# Groups A1 x B1 #
##################

# Permutation test.
perm_t_GA1xB1_un <- perm_test(hm_dist_un_m[,"A1"], hm_dist_un_m[,"B1"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1xB1_un$pval

# Observed correlation.
perm_t_GA1xB1_un$obs_corr

# Scatter plot.
plot(hm_dist_un_m[,"A1"], hm_dist_un_m[,"B1"], xlab = "A1", ylab = "B1", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1xB1_un$perm_corr, main = "A1 x B1", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1xB1_un$obs_corr)*1.3),2),round(((perm_t_GA1xB1_un$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA1xB1_un$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups A1 x B2 #
##################

# Permutation test.
perm_t_GA1xB2_un <- perm_test(hm_dist_un_m[,"A1"], hm_dist_un_m[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1xB2_un$pval

# Observed correlation.
perm_t_GA1xB2_un$obs_corr

# Scatter plot.
plot(hm_dist_un_m[,"A1"], hm_dist_un_m[,"B2"], xlab = "A1", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1xB2_un$perm_corr, main = "A1 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1xB2_un$obs_corr)*1.3),2),round(((perm_t_GA1xB2_un$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA1xB2_un$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups A2 x B1 #
##################

# Permutation test.
perm_t_GA2xB1_un <- perm_test(hm_dist_un_m[,"A2"], hm_dist_un_m[,"B1"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA2xB1_un$pval

# Observed correlation.
perm_t_GA2xB1_un$obs_corr

# Scatter plot.
plot(hm_dist_un_m[,"A2"], hm_dist_un_m[,"B1"], xlab = "A2", ylab = "B1", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA2xB1_un$perm_corr, main = "A2 x B1", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA2xB1_un$obs_corr)*1.3),2),round(((perm_t_GA2xB1_un$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA2xB1_un$obs_corr, lty = 2, col = "red")
dev.off()

#----------#

##################
# Groups A2 x B2 #
##################

# Permutation test.
perm_t_GA2xB2_un <- perm_test(hm_dist_un_m[,"A2"], hm_dist_un_m[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA2xB2_un$pval

# Observed correlation.
perm_t_GA2xB2_un$obs_corr

# Scatter plot.
plot(hm_dist_un_m[,"A2"], hm_dist_un_m[,"B2"], xlab = "A2", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA2xB2_un$perm_corr, main = "A2 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA2xB2_un$obs_corr)*1.3),2),round(((perm_t_GA2xB2_un$obs_corr)*1.3),3)), breaks = 50)
abline(v = perm_t_GA2xB2_un$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------------------------------------#

# Save data.
save(perm_t_GA1x2_un, perm_t_GB1x2_un, perm_t_GA1xB1_un, perm_t_GA1xB2_un, perm_t_GA2xB1_un, perm_t_GA2xB2_un,
     file = "data_out/permut_tests_uncond.RDA")

#--------------------------------------------------#

#########################
### MORLON'S APPROACH ###
#########################

# Set list to re-scaled trees.
tree_set_un <- setNames(vector(mode = "list", length = dim(hm_dist_un_m)[2]), colnames(hm_dist_un_m))
tree_set_un

j = 1
# Loop over.
for (j in 1:length(tree_set_un)) {
  
  # Copy original tree.
  tree_c <- tree
  
  # Re-scale branches by mean observed changes. Replace branch lengths from original tree.
  tree_c$edge.length <- hm_dist_un_m[,j]
  #tree_c$edge.length[tree_c$edge.length == 0] <- 1e-3
  
  # Store re-scaled tree.
  tree_set_un[[j]] <- tree_c
  
}

# Load packages.
library(RPANDA)

# Calculate JSD among trees. Probably a bad idea with too few trees!
jsd_un = JSDtree(tree_set_un, meth = c("standard"))

# Get clusters.
JSDtree_cluster(jsd_un)

# Notes:
# Probably not working properly since several branches have 0 length.
