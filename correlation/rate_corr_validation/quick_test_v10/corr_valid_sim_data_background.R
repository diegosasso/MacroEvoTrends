
# Load libraries.
library(ape)
library(phytools)
library(rphenoscate)
library(ontophylo)
library(corHMM)

# Import some functions.
source("R/utils.R")
source("R/helpers.R")

# Load data.
load("data_out/sim_traits.RDA")

#--------------------------------------------------#

#######################
# Reorganize datasets #
#######################

# OPTION 1 # 

# Remove only noisy traits.
stm = list("A1" = stm$A1[c(1:15)],
           "A2" = stm$A2[c(1:15)],
           "B1" = stm$B1[c(1:15)],
           "B2" = stm$B2[c(1:15)])

# Save data.
saveRDS(stm, "data_out/stm_alt1.RDS")

# OPTION 2 # 

# Remove traits with shift regimes on focal branches (i.e., keep only noisy traits).
stm = list("A1" = stm$A1[c(16:30)],
           "A2" = stm$A2[c(16:30)],
           "B1" = stm$B1[c(16:30)],
           "B2" = stm$B2[c(16:30)])

# Save data.
saveRDS(stm, "data_out/stm_alt2.RDS")

#--------------------------------------------------#

##############################
# Amalgamate stochastic maps #
##############################

# Set parameters.
nstm = 200
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
stm_amalg_alt <- list("A1" = stm_GA1_amalg,
                      "A2" = stm_GA2_amalg,
                      "B1" = stm_GB1_amalg,
                      "B2" = stm_GB2_amalg)

# Save data.
saveRDS(stm_amalg_alt, "data_out/stm_amalg_alt1.RDS")
saveRDS(stm_amalg_alt, "data_out/stm_amalg_alt2.RDS")

#--------------------------------------------------#

####################
# Get branch rates #
####################

# Merge tree bins.
stm_merg <- lapply(stm_amalg_alt, function(x) merge_tree_cat_list(x) )
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
     border = F, xlim = c(-round(((perm_t_GA1xB1$obs_corr)*5.0),2),round(((perm_t_GA1xB1$obs_corr)*5.0),3)), breaks = 50)
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
     border = F, xlim = c(-round(((perm_t_GA1xB2$obs_corr)*15.0),2),round(((perm_t_GA1xB2$obs_corr)*15.0),3)), breaks = 50)
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
     border = F, xlim = c(-round(((perm_t_GA2xB1$obs_corr)*15.0),2),round(((perm_t_GA2xB1$obs_corr)*15.0),3)), breaks = 50)
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
     border = F, xlim = c(-round(((perm_t_GA2xB2$obs_corr)*10.0),2),round(((perm_t_GA2xB2$obs_corr)*10.0),3)), breaks = 50)
abline(v = perm_t_GA2xB2$obs_corr, lty = 2, col = "red")
dev.off()

#--------------------------------------------------#

# NOTES:
#
# When comparing these results with the ones obtained from the full analyses
# (i.e., traits with high rates on focal branches + noisy traits)
# it seems that we have lower p-values and considerably higher correlations
# between groups A1 x B1, A1 x B2, A2 x B1, and A2 x B2, which should not be correlated at all;
# I suspect the mere fact of having some branches with very high rates when compared to the background,
# irrespective if they are the same or different branches, is enough to increase the chances of false positives;
# probably the problem gets worse when more traits are amalgamated;
# Also, only "good" traits or only "noisy" traits produce good p-values
# (i.e., significant for A1 x A2 and B1 x B2 but not otherwise);
# When "good" + "noisy" are amalgamated together (groups with 30 traits total),
# p-values for non-correlated partitions decreases significantly, almost reaching significance levels (false positives!)
# Maybe this is just a matter of adding more traits (15 vs. 30);
# Should we re-scale rates using mean changes per character?



