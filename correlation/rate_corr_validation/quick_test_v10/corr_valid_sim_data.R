
#--------------------------------------------------------------------------------#

# EXPLANATION NOTES:
#
# This code is used to simulate two groups of binary traits.
#
# Traits will be simulated under a two-regime model:
# Regime 1: a background regime based on a symmetric Q matrix (Q_base);
# Regime 2: a shift regime based on a highly asymmetric Q matrix (Q_jump);
#
# Regime 1 applies to the whole tree (except focal branch);
# Regime 2 applies only to a focal branch;
#
# The shift regime is intended to produce irreversible transitions on focal branches of interest;
# The irreversible transitions result in accumulation of directional changes on focal branches when several traits are amalgamated;
# Such changes mimic "synapomorphies" as seem in typical phylogenetic datasets and result in higher evolutionary rates on the focal branches;
#
# The first group (A) comprises traits simulated to produce high evolutionary rates on branches 4, 5, and 6 of the Hymenoptera tree;
# The first group (B) comprises traits simulated to produce high evolutionary rates on 3 other random branches (i.e., excluding 4-6;
# To produce simulated datasets similar to typical phylogenetic datasets, branches leading to tips were removed from the sampling pool
# since autapomorphies are usually note included in empirical datasets, only parsimony-informative traits (i.e., variable in two or more taxa);
#
# For each focal branch, 5 binary traits were simulated thus resulting in 15 traits in total for each group;
# Since not all traits in empirical phylogenetic datasets produce unique state-transitions, for each group,
# 15 additional binary traits were simulated under the background regime only and added to the initial traits
# thus resulting in 30 binary traits per group.
#
# For each simulated trait in each group, models were fitted using corHMM and then 200 stochastic maps were sampled;
# Then, all 30 traits in each group were amalgamated using the PARAMO pipeline 
# resulting in 200 amalgamated stochastic maps for each group;
#
# Finally, groups A and B were split in half, A1 and A2, B1 and B2, respectively, with 100 stochastic maps each;
# Since in each group (A and B) all traits were simulated under the same generative model (Q_base/Q_jump),
# same tip state distribution, same fitted models from corHMM, and same focal branches,
# just differing in the sampled stochastic histories in each sample of 100 maps,
# we have the following expectations:
# (1) branch rates in A1 x A2 will be correlated;
# (2) branch rates in B1 x B2 will be correlated;
# (3) branch rates in A1 x B1, A1 x B2, A2 x B1, and A2 x B2 will not be correlated;
#
# The core goal of this script is to answer the question:
# Can we detect TRUE BRANCH RATE CORRELATIONS using common phylogenetic datasets when only a few branches have signal for correlated rates?

#--------------------------------------------------------------------------------#

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

#--------------------------------------------------#

##############
# Set models #
##############

# Set base Q-matrix.
Q_base <- initQ(c(0, 1), c(0.001,0.001))
Q_base

# Set Q-matrix of shift.
Q_jump <- initQ(c(0, 1), c(10, 1e-10))
Q_jump

#--------------------------------------------------#

# Set some parameters.
nsim = 5
nstm = 200

#--------------------------------------------------#

#################
# Simulate data #
#################

# Define a set of focal branches for the shifts ("synapomorphies").
br = which(!hym_tree$edge[,2] %in% c(1:length(hym_tree$tip.label)))

# Note: this set excludes all terminal branches since phylo data
# usually does not include autapomorphies.

# Define two sets of branches (3 branches each) to simulate correlated rates.
br_s1 <- c(4,5,6) # series of branches leading to Apocrita.
br_s2 <- sample(br[-c(4,5,6)],3) # random sample of other branches.

br_sample = c(br_s1, br_s2)

#--------------------#

# Simulate tip data.

# Create folders to store data.
suppressWarnings(dir.create("data_out"))

# Set a list to store simulated tip data.
sim_data <- setNames(vector(mode = "list", length = length(br_sample)), paste0("BR_", br_sample))

# Simulate phylogenetically-informative tip data.
for (i in 1:length(br_sample)) {
  
  # Simulate tip data.
  sim_data[[i]] <- sim_info_data(tree = hym_tree, focal_edge = br_sample[[i]], Q_base = Q_base, Q_jump = Q_jump, nsim = nsim, jump = T)
  
}

# Set a list to store simulated noisy tip data.
sim_data_n <- vector(mode = "list", length = length(br_sample))

# Simulate noisy tip data.
for (i in 1:length(br_sample)) {
  
  # Simulate tip data.
  sim_data_n[[i]] <- sim_info_data(tree = hym_tree, Q_base = Q_base, nsim = nsim, jump = F)
  
}

#--------------------#

# Quick checks.
cbind(sim_data[[1]][,1], sim_data_n[[1]][,1])
cbind(sim_data[[2]][,1], sim_data_n[[2]][,1])
cbind(sim_data[[3]][,1], sim_data_n[[3]][,1])

#--------------------#

# Set a list to store sampled character histories.
stm_data <- setNames(vector(mode = "list", length = length(br_sample)), paste0("BR_", br_sample))

# Simulate data and sample character histories.
for (i in 1:length(br_sample)) {
  
  # Sample character histories conditioned on tip data.
  stm_data[[i]] <- sample_cond(hym_tree, sim_data[[i]], nsim, nstm)
  
}

# Set a list to store sampled character histories (noisy data).
stm_data_n <- vector(mode = "list", length = length(br_sample))

# Simulate data and sample character histories.
for (i in 1:length(br_sample)) {
  
  # Sample character histories conditioned on tip data.
  stm_data_n[[i]] <- sample_cond(hym_tree, sim_data_n[[i]], nsim, nstm)
  
}

#--------------------#

# Quick checks.
densityMap(stm_data$BR_4$smaps$C1[1:100], fsize = 0.8)
densityMap(stm_data$BR_4$smaps$C2[1:100], fsize = 0.8)
densityMap(stm_data$BR_5$smaps$C1[1:100], fsize = 0.8)
densityMap(stm_data$BR_5$smaps$C2[1:100], fsize = 0.8)
densityMap(stm_data$BR_6$smaps$C1[1:100], fsize = 0.8)
densityMap(stm_data$BR_6$smaps$C2[1:100], fsize = 0.8)

densityMap(stm_data$BR_12$smaps$C1[1:100], fsize = 0.8)
densityMap(stm_data$BR_12$smaps$C2[1:100], fsize = 0.8)
densityMap(stm_data$BR_78$smaps$C1[1:100], fsize = 0.8)
densityMap(stm_data$BR_78$smaps$C2[1:100], fsize = 0.8)
densityMap(stm_data$BR_15$smaps$C1[1:100], fsize = 0.8)
densityMap(stm_data$BR_15$smaps$C2[1:100], fsize = 0.8)

densityMap(stm_data_n[[1]]$smaps$C1[1:100], fsize = 0.8)
densityMap(stm_data_n[[1]]$smaps$C2[1:100], fsize = 0.8)
densityMap(stm_data_n[[1]]$smaps$C1[1:100], fsize = 0.8)

#--------------------------------------------------#

#################################
# Merge and reorganize datasets #
#################################

# Organize simulated tip data.
sim_traits <- list("A" = do.call(cbind, c(sim_data[c(1:3)], sim_data_n[c(1:3)])),
                   "B" = do.call(cbind, c(sim_data[c(4:6)], sim_data_n[c(4:6)])))

colnames(sim_traits$A) <- paste0("C", 1:length(colnames(sim_traits$A)))
colnames(sim_traits$B) <- paste0("C", 1:length(colnames(sim_traits$B)))
sim_traits

#--------------------#

# Organize fitted models maps.
fitted_models <- list("A" = c(stm_data$BR_4$models,
                               stm_data$BR_5$models,
                               stm_data$BR_6$models,
                               stm_data_n[[1]]$models,
                               stm_data_n[[2]]$models,
                               stm_data_n[[3]]$models),

                      "B" = c(stm_data$BR_12$models,
                               stm_data$BR_78$models,
                               stm_data$BR_15$models,
                               stm_data_n[[4]]$models,
                               stm_data_n[[5]]$models,
                               stm_data_n[[6]]$models)
                      )

fitted_models <- lapply(fitted_models, function(x) setNames(x, paste0("C", 1:length(x))) )
fitted_models

#--------------------#

# Organize stochastic maps.
stm <- list("A1" = c(lapply(stm_data$BR_4$smaps, function(x) x[1:100]),
                     lapply(stm_data$BR_5$smaps, function(x) x[1:100]),
                     lapply(stm_data$BR_6$smaps, function(x) x[1:100]),
                     lapply(stm_data_n[[1]]$smaps, function(x) x[1:100]),
                     lapply(stm_data_n[[2]]$smaps, function(x) x[1:100]),
                     lapply(stm_data_n[[3]]$smaps, function(x) x[1:100])),
            
            "A2" = c(lapply(stm_data$BR_4$smaps, function(x) x[101:200]),
                     lapply(stm_data$BR_5$smaps, function(x) x[101:200]),
                     lapply(stm_data$BR_6$smaps, function(x) x[101:200]),
                     lapply(stm_data_n[[1]]$smaps, function(x) x[101:200]),
                     lapply(stm_data_n[[2]]$smaps, function(x) x[101:200]),
                     lapply(stm_data_n[[3]]$smaps, function(x) x[101:200])),
            
            "B1" = c(lapply(stm_data$BR_12$smaps, function(x) x[1:100]),
                     lapply(stm_data$BR_78$smaps, function(x) x[1:100]),
                     lapply(stm_data$BR_15$smaps, function(x) x[1:100]),
                     lapply(stm_data_n[[4]]$smaps, function(x) x[1:100]),
                     lapply(stm_data_n[[5]]$smaps, function(x) x[1:100]),
                     lapply(stm_data_n[[6]]$smaps, function(x) x[1:100])),
            
            "B2" = c(lapply(stm_data$BR_12$smaps, function(x) x[101:200]),
                     lapply(stm_data$BR_78$smaps, function(x) x[101:200]),
                     lapply(stm_data$BR_15$smaps, function(x) x[101:200]),
                     lapply(stm_data_n[[4]]$smaps, function(x) x[101:200]),
                     lapply(stm_data_n[[5]]$smaps, function(x) x[101:200]),
                     lapply(stm_data_n[[6]]$smaps, function(x) x[101:200]))
            )

stm <- lapply(stm, function(x) setNames(x, paste0("C", 1:length(x))) )
stm

#--------------------#

# Save data.
save(hym_tree, br_s1, br_s2, Q_base, Q_jump, sim_traits, fitted_models, stm, file = "data_out/sim_traits.RDA")


#--------------------------------------------------#

##############################
# Amalgamate stochastic maps #
##############################

# Set parameters.
res = 200

# Discretize maps.
stm_GA1_discr <- lapply(stm$A1, function(x) discr_Simmap_all(x, res) )
stm_GA2_discr <- lapply(stm$A2, function(x) discr_Simmap_all(x, res) )
stm_GB1_discr <- lapply(stm$B1, function(x) discr_Simmap_all(x, res) )
stm_GB2_discr <- lapply(stm$B2, function(x) discr_Simmap_all(x, res) )

# Amalgamate maps.
stm_GA1_amalg <- paramo.list(names(stm_GA1_discr), tree.list = stm_GA1_discr, ntrees = (nstm/2) )
stm_GA2_amalg <- paramo.list(names(stm_GA2_discr), tree.list = stm_GA2_discr, ntrees = (nstm/2) )
stm_GB1_amalg <- paramo.list(names(stm_GB1_discr), tree.list = stm_GB1_discr, ntrees = (nstm/2) )
stm_GB2_amalg <- paramo.list(names(stm_GB2_discr), tree.list = stm_GB2_discr, ntrees = (nstm/2) )

# Merge datasets.
stm_amalg <- list("A1" = stm_GA1_amalg,
                  "A2" = stm_GA2_amalg,
                  "B1" = stm_GB1_amalg,
                  "B2" = stm_GB2_amalg)

# Save data.
saveRDS(stm_amalg, "data_out/stm_amalg.RDS")

#--------------------------------------------------#

####################
# Get branch rates #
####################

# Import data.
stm_amalg <- readRDS("data_out/stm_amalg.RDS")

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

#--------------------------------------------------#

########################
# Possible corrections #
########################

# Try: Removing branches leading to tips.
br_out <- which(!hym_tree$edge[,2] %in% c(1:length(hym_tree$tip.label)))
#hm_dist_m_rt <- hm_dist_m_rt[br_out,]
#hm_dist_m_rt 

# Try: residuals from regression between number of changes and branch lengths.
res_regr <- cbind(lm(hm_dist_m[,1] ~ hym_tree$edge.length)$residuals,
                  lm(hm_dist_m[,2] ~ hym_tree$edge.length)$residuals,
                  lm(hm_dist_m[,3] ~ hym_tree$edge.length)$residuals,
                  lm(hm_dist_m[,4] ~ hym_tree$edge.length)$residuals)
colnames(res_regr) <- colnames(hm_dist_m_rt)
hm_dist_m_rt <- res_regr
hm_dist_m_rt <- hm_dist_m_rt[br_out,]
hm_dist_m_rt 

#--------------------------------------------------#

#####################
# PERMUTATION TESTS #
#####################

#################### RAW RATES ####################

########################
# Should be correlated #
########################

##################
# Groups A1 x A2 #
##################

# Permutation test.
perm_t_GA1x2 <- perm_test(hm_dist_m_rt[,"A1"], hm_dist_m_rt[,"A2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1x2$pval

# Observed correlation.
perm_t_GA1x2$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"A1"], hm_dist_m_rt[,"A2"], xlab = "A1", ylab = "A2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1x2$perm_corr, main = "A1 x A2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1x2$obs_corr)*1.3),2),round(((perm_t_GA1x2$obs_corr)*1.3),2)), breaks = 50)
abline(v = perm_t_GA1x2$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------#

##################
# Groups B1 x B2 #
##################

# Permutation test.
perm_t_GB1x2 <- perm_test(hm_dist_m_rt[,"B1"], hm_dist_m_rt[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GB1x2$pval

# Observed correlation.
perm_t_GB1x2$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"B1"], hm_dist_m_rt[,"B2"], xlab = "B1", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GB1x2$perm_corr, main = "B1 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GB1x2$obs_corr)*1.3),2),round(((perm_t_GB1x2$obs_corr)*1.3),2)), breaks = 50)
abline(v = perm_t_GB1x2$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------------------------------------#

############################
# Should not be correlated #
############################

##################
# Groups A1 x B1 #
##################

# Permutation test.
perm_t_GA1xB1 <- perm_test(hm_dist_m_rt[,"A1"], hm_dist_m_rt[,"B1"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1xB1$pval

# Observed correlation.
perm_t_GA1xB1$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"A1"], hm_dist_m_rt[,"B1"], xlab = "A1", ylab = "B1", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1xB1$perm_corr, main = "A1 x B1", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1xB1$obs_corr)*3.0),2),round(((perm_t_GA1xB1$obs_corr)*3.0),3)), breaks = 50)
abline(v = perm_t_GA1xB1$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------#

##################
# Groups A1 x B2 #
##################

# Permutation test.
perm_t_GA1xB2 <- perm_test(hm_dist_m_rt[,"A1"], hm_dist_m_rt[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA1xB2$pval

# Observed correlation.
perm_t_GA1xB2$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"A1"], hm_dist_m_rt[,"B2"], xlab = "A1", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA1xB2$perm_corr, main = "A1 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA1xB2$obs_corr)*3.0),2),round(((perm_t_GA1xB2$obs_corr)*3.0),3)), breaks = 50)
abline(v = perm_t_GA1xB2$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------#

##################
# Groups A2 x B1 #
##################

# Permutation test.
perm_t_GA2xB1 <- perm_test(hm_dist_m_rt[,"A2"], hm_dist_m_rt[,"B1"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA2xB1$pval

# Observed correlation.
perm_t_GA2xB1$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"A2"], hm_dist_m_rt[,"B1"], xlab = "A2", ylab = "B1", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA2xB1$perm_corr, main = "A2 x B1", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA2xB1$obs_corr)*3.0),2),round(((perm_t_GA2xB1$obs_corr)*3.0),3)), breaks = 50)
abline(v = perm_t_GA2xB1$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------#

##################
# Groups A2 x B2 #
##################

# Permutation test.
perm_t_GA2xB2 <- perm_test(hm_dist_m_rt[,"A2"], hm_dist_m_rt[,"B2"], nperm = 10000, method = "spearman")

# P-value.
perm_t_GA2xB2$pval

# Observed correlation.
perm_t_GA2xB2$obs_corr

# Scatter plot.
plot(hm_dist_m_rt[,"A2"], hm_dist_m_rt[,"B2"], xlab = "A2", ylab = "B2", pch = 19, col = "red")
dev.off()

# Histogram.
hist(perm_t_GA2xB2$perm_corr, main = "A2 x B2", xlab = "correlation", 
     border = F, xlim = c(-round(((perm_t_GA2xB2$obs_corr)*3.0),2),round(((perm_t_GA2xB2$obs_corr)*3.0),3)), breaks = 50)
abline(v = perm_t_GA2xB2$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------------------------------------#

