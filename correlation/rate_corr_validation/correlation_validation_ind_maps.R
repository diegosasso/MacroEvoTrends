
# Load libraries.
library(ape)
library(phytools)

# Import some functions.
source("R/utils.R")
source("R/helpers.R")

################################
### TEST FOR INDIVIDUAL MAPS ###
################################

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

# Convert to rates.
hm_dist <- lapply(hm_dist, function(x) lapply(x, function(y) y/tree$edge.length ) )

# Set lists to store results.
GA1x2_corr_ind <- numeric()
GA1x2_pval_ind <- numeric()
GB1x2_corr_ind <- numeric()
GB1x2_pval_ind <- numeric()
GA1xB1_corr_ind <- numeric()
GA1xB1_pval_ind <- numeric()
GA1xB2_corr_ind <- numeric()
GA1xB2_pval_ind <- numeric()
GA2xB1_corr_ind <- numeric()
GA2xB1_pval_ind <- numeric()
GA2xB2_corr_ind <- numeric()
GA2xB2_pval_ind <- numeric()

# Loop over all stochastic maps.

#----------#

## A1 x A2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist$A1[[i]], hm_dist$A2[[i]], nperm = 10000, method = "spearman")
  GA1x2_corr_ind[[i]] <- perm_t$obs_corr
  GA1x2_pval_ind[[i]] <- perm_t$pval  
  
}

#----------#

## B1 x B2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist$B1[[i]], hm_dist$B2[[i]], nperm = 10000, method = "spearman")
  GB1x2_corr_ind[[i]] <- perm_t$obs_corr
  GB1x2_pval_ind[[i]] <- perm_t$pval
  
}

#----------#

## A1 x B1 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist$A1[[i]], hm_dist$B1[[i]], nperm = 10000, method = "spearman")
  GA1xB1_corr_ind[[i]] <- perm_t$obs_corr
  GA1xB1_pval_ind[[i]] <- perm_t$pval
  
}

#----------#

## A1 x B2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist$A1[[i]], hm_dist$B2[[i]], nperm = 10000, method = "spearman")
  GA1xB2_corr_ind[[i]] <- perm_t$obs_corr
  GA1xB2_pval_ind[[i]] <- perm_t$pval
  
}

#----------#

## A2 x B1 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist$A2[[i]], hm_dist$B1[[i]], nperm = 10000, method = "spearman")
  GA2xB1_corr_ind[[i]] <- perm_t$obs_corr
  GA2xB1_pval_ind[[i]] <- perm_t$pval
  
  
}

#----------#

## A2 x B2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist$A2[[i]], hm_dist$B2[[i]], nperm = 10000, method = "spearman")
  GA2xB2_corr_ind[[i]] <- perm_t$obs_corr
  GA2xB2_pval_ind[[i]] <- perm_t$pval  
  
}

#----------#

# Check proportion of significant p-values.
sum(GA1x2_pval_ind <= 0.05)/100
sum(GB1x2_pval_ind <= 0.05)/100

sum(GA1xB1_pval_ind <= 0.05)/100
sum(GA1xB2_pval_ind <= 0.05)/100

sum(GA2xB1_pval_ind <= 0.05)/100
sum(GA2xB2_pval_ind <= 0.05)/100

#----------#

#--------------------------------------------------#



################################################################################



#--------------------------------------------------#
# UNCONDITIONED ON DATA
#--------------------------------------------------#

# Import data.
stm_un_amalg <- readRDS("data_out/stm_amalg_uncond.RDS")

# Merge tree bins.
stm_un_merg <- lapply(stm_un_amalg, function(x) merge_tree_cat_list(x) )

# Get max Hamming distances.
hm_dist_un <- lapply(stm_un_merg, function(x) hm_dist_br(x) )

# Convert to rates.
hm_dist_un <- lapply(hm_dist_un, function(x) lapply(x, function(y) y/tree$edge.length ) )

# Set lists to store results.
GA1x2_un_corr_ind <- numeric()
GA1x2_un_pval_ind <- numeric()
GB1x2_un_corr_ind <- numeric()
GB1x2_un_pval_ind <- numeric()
GA1xB1_un_corr_ind <- numeric()
GA1xB1_un_pval_ind <- numeric()
GA1xB2_un_corr_ind <- numeric()
GA1xB2_un_pval_ind <- numeric()
GA2xB1_un_corr_ind <- numeric()
GA2xB1_un_pval_ind <- numeric()
GA2xB2_un_corr_ind <- numeric()
GA2xB2_un_pval_ind <- numeric()

# Loop over all stochastic maps.

#----------#

## A1 x A2 ##

for (i in 1:100) {

  perm_t <- perm_test(hm_dist_un$A1[[i]], hm_dist_un$A2[[i]], nperm = 10000, method = "spearman")
  GA1x2_un_corr_ind[[i]] <- perm_t$obs_corr
  GA1x2_un_pval_ind[[i]] <- perm_t$pval  
  
}

#----------#

## B1 x B2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_un$B1[[i]], hm_dist_un$B2[[i]], nperm = 10000, method = "spearman")
  GB1x2_un_corr_ind[[i]] <- perm_t$obs_corr
  GB1x2_un_pval_ind[[i]] <- perm_t$pval
  
}

#----------#

## A1 x B1 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_un$A1[[i]], hm_dist_un$B1[[i]], nperm = 10000, method = "spearman")
  GA1xB1_un_corr_ind[[i]] <- perm_t$obs_corr
  GA1xB1_un_pval_ind[[i]] <- perm_t$pval
  
}

#----------#

## A1 x B2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_un$A1[[i]], hm_dist_un$B2[[i]], nperm = 10000, method = "spearman")
  GA1xB2_un_corr_ind[[i]] <- perm_t$obs_corr
  GA1xB2_un_pval_ind[[i]] <- perm_t$pval
  
}

#----------#

## A2 x B1 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_un$A2[[i]], hm_dist_un$B1[[i]], nperm = 10000, method = "spearman")
  GA2xB1_un_corr_ind[[i]] <- perm_t$obs_corr
  GA2xB1_un_pval_ind[[i]] <- perm_t$pval
  
  
}

#----------#

## A2 x B2 ##

for (i in 1:100) {

  perm_t <- perm_test(hm_dist_un$A2[[i]], hm_dist_un$B2[[i]], nperm = 10000, method = "spearman")
  GA2xB2_un_corr_ind[[i]] <- perm_t$obs_corr
  GA2xB2_un_pval_ind[[i]] <- perm_t$pval  
  
}

#----------#

# Check proportion of significant p-values.
sum(GA1x2_un_pval_ind <= 0.05)/100
sum(GB1x2_un_pval_ind <= 0.05)/100

sum(GA1xB1_un_pval_ind <= 0.05)/100
sum(GA1xB2_un_pval_ind <= 0.05)/100

sum(GA2xB1_un_pval_ind <= 0.05)/100
sum(GA2xB2_un_pval_ind <= 0.05)/100

#----------#

#--------------------------------------------------#
# DIFFERENCE DATA
#--------------------------------------------------#

# Get only difference in changes
#hm_dist_c <- list("A1" = mapply(x = hm_dist$A1, y = hm_dist_un$A1, function(x,y) sqrt((x - y)^2), SIMPLIFY = F),
#                  "A2" = mapply(x = hm_dist$A2, y = hm_dist_un$A2, function(x,y) sqrt((x - y)^2), SIMPLIFY = F),
#                  "B1" = mapply(x = hm_dist$B1, y = hm_dist_un$B1, function(x,y) sqrt((x - y)^2), SIMPLIFY = F),
#                  "B2" = mapply(x = hm_dist$B2, y = hm_dist_un$B2, function(x,y) sqrt((x - y)^2), SIMPLIFY = F))

hm_dist_c <- list("A1" = mapply(x = hm_dist$A1, y = hm_dist_un$A1, function(x,y) x-y, SIMPLIFY = F),
                  "A2" = mapply(x = hm_dist$A2, y = hm_dist_un$A2, function(x,y) x-y, SIMPLIFY = F),
                  "B1" = mapply(x = hm_dist$B1, y = hm_dist_un$B1, function(x,y) x-y, SIMPLIFY = F),
                  "B2" = mapply(x = hm_dist$B2, y = hm_dist_un$B2, function(x,y) x-y, SIMPLIFY = F))

# Set lists to store results.
GA1x2_c_corr_ind <- numeric()
GA1x2_c_pval_ind <- numeric()
GB1x2_c_corr_ind <- numeric()
GB1x2_c_pval_ind <- numeric()
GA1xB1_c_corr_ind <- numeric()
GA1xB1_c_pval_ind <- numeric()
GA1xB2_c_corr_ind <- numeric()
GA1xB2_c_pval_ind <- numeric()
GA2xB1_c_corr_ind <- numeric()
GA2xB1_c_pval_ind <- numeric()
GA2xB2_c_corr_ind <- numeric()
GA2xB2_c_pval_ind <- numeric()

# Loop over all stochastic maps.

#----------#

## A1 x A2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_c$A1[[i]], hm_dist_c$A2[[i]], nperm = 10000, method = "spearman")
  GA1x2_c_corr_ind[[i]] <- perm_t$obs_corr
  GA1x2_c_pval_ind[[i]] <- perm_t$pval  
  
}

#----------#

## B1 x B2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_c$B1[[i]], hm_dist_c$B2[[i]], nperm = 10000, method = "spearman")
  GB1x2_c_corr_ind[[i]] <- perm_t$obs_corr
  GB1x2_c_pval_ind[[i]] <- perm_t$pval
  
}

#----------#

## A1 x B1 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_c$A1[[i]], hm_dist_c$B1[[i]], nperm = 10000, method = "spearman")
  GA1xB1_c_corr_ind[[i]] <- perm_t$obs_corr
  GA1xB1_c_pval_ind[[i]] <- perm_t$pval
  
}

#----------#

## A1 x B2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_c$A1[[i]], hm_dist_c$B2[[i]], nperm = 10000, method = "spearman")
  GA1xB2_c_corr_ind[[i]] <- perm_t$obs_corr
  GA1xB2_c_pval_ind[[i]] <- perm_t$pval
  
}

#----------#

## A2 x B1 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_c$A2[[i]], hm_dist_c$B1[[i]], nperm = 10000, method = "spearman")
  GA2xB1_c_corr_ind[[i]] <- perm_t$obs_corr
  GA2xB1_c_pval_ind[[i]] <- perm_t$pval
  
  
}

#----------#

## A2 x B2 ##

for (i in 1:100) {
  
  perm_t <- perm_test(hm_dist_c$A2[[i]], hm_dist_c$B2[[i]], nperm = 10000, method = "spearman")
  GA2xB2_c_corr_ind[[i]] <- perm_t$obs_corr
  GA2xB2_c_pval_ind[[i]] <- perm_t$pval  
  
}

#----------#

# Check proportion of significant p-values.
sum(GA1x2_c_pval_ind <= 0.05)/100
sum(GB1x2_c_pval_ind <= 0.05)/100

sum(GA1xB1_c_pval_ind <= 0.05)/100
sum(GA1xB2_c_pval_ind <= 0.05)/100

sum(GA2xB1_c_pval_ind <= 0.05)/100
sum(GA2xB2_c_pval_ind <= 0.05)/100

#----------#

#--------------------------------------------------#
