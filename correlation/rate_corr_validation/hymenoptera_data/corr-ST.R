# Load libraries.
library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)

# Import some functions.
source("R/helpers.R")
source("R/utils-ST.R")

# Load Hymenoptera dataset.
#load("data_hym/paramo_stm_final.RDA")

  
#########################
# Check empirical rates #
#########################

# Import data.
stm_amalg <- readRDS("data_hym/anat_ent.RDS")
head(stm_amalg)
length(stm_amalg$mesopectus)

# Merge tree bins.
stm_merg <- lapply(stm_amalg, function(x) merge_tree_cat_list(x) )
#saveRDS(stm_merg, "data_hym/stm_merg.RDS")
#get_branch_rate_sample_smaps(stm_merg$cranium)
#get_branch_rate_across_smaps(stm_merg$cranium)

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

#------- Calculate pairwise correlation 

stm_merg <- readRDS("data_hym/stm_merg.RDS")

result <-test_correlation_with_resampling(
    stm_merg$mesopectus, stm_merg$`female genitalia`,
    n.samples = 1000,
    cor_method = "spearman",
    remove_zeros = FALSE,
    scale = FALSE,
    log.trans = TRUE,
    log.constant = 1e-6,
    quantile.trim = c(0.05, 0.95)
)

print(result$pval)
print(result$obs_corr)
hist(result$null_distribution)

names(stm_merg)
cat(names(stm_merg), sep='\n')

corr.spear <- run_all_pairwise_correlations(
  stm_merg[1:4],
  n.samples = 1000,
  cor_method = "spearman",
  remove_zeros = FALSE,
  scale = FALSE,
  log.trans = F,
  log.constant = 1e-6,
  quantile.trim = c(0.05, 0.95)
)


corr.spear <-run_all_pairwise_correlations_parallel(
  stm_merg,
  n.samples = 1000,
  cor_method = "spearman",
  remove_zeros = FALSE,
  scale = FALSE,
  log.trans = F,
  log.constant = 1e-6,
  quantile.trim = c(0.01, 0.99),
  n.cores = parallel::detectCores() - 1
)

corr.spear.noZ <-run_all_pairwise_correlations_parallel(
  stm_merg,
  n.samples = 1000,
  cor_method = "spearman",
  remove_zeros = T,
  scale = FALSE,
  log.trans = F,
  log.constant = 1e-6,
  quantile.trim = c(0.01, 0.99),
  n.cores = parallel::detectCores() - 1
)

corr.pearson.log <-run_all_pairwise_correlations_parallel(
  stm_merg,
  n.samples = 1000,
  cor_method = "pearson",
  remove_zeros = FALSE,
  scale = FALSE,
  log.trans = T,
  log.constant = 1e-6,
  quantile.trim = c(0.01, 0.99),
  n.cores = parallel::detectCores() - 1
)

corr.pearson.scaled <-run_all_pairwise_correlations_parallel(
  stm_merg,
  n.samples = 1000,
  cor_method = "pearson",
  remove_zeros = FALSE,
  scale = T,
  log.trans = F,
  log.constant = 1e-6,
  quantile.trim = c(0.01, 0.99),
  n.cores = parallel::detectCores() - 1
)



pval(corr.spear)
pval(corr.spear, use.bonf = T, 0.05)
plot(pval(corr.spear), pval(corr.spear.noZ))
plot(pval(corr.spear), pval(corr.spear.noZ), xlim = c(0, 0.05))
plot(pval(corr.spear), pval(corr.pearson.log), xlim = c(0, 0.05))


#--- corellation, rough, not filtered
#pv <- pval(corr.pearson.scaled, use.bonf = T) 
pv <- pval(corr.pearson.log, use.bonf = T, 0.05) 
cv <- corr_val(corr.pearson.log)
hist(cv)
cv[names(pv)]


#cv <- corr_val(corr.spear)

# branch rate
rates_list <- lapply(stm_merg, get_branch_rate_across_smaps)
rates <- do.call(cbind, rates_list)

head(rates)
rate_total <- apply(rates, 2, sum)

# Create summed rates for each pair in `cv`
rate_sums <- sapply(names(cv), function(pair) {
  parts <- strsplit(pair, " vs ")[[1]]
  sum(rate_total[parts])
  #mean(rate_total[parts])
})

plot(cv, rate_sums)
cor.test(rate_sums, cv, method = "spearman")

plot(rate_sums, cv,
     xlab = "Sum of Total Rates",
     ylab = "Pairwise Correlation (cv)",
     main = "Correlation of Total Rates vs Pairwise Correlations")
abline(lm(cv ~ rate_sums), col = "blue")


sig_idx <- pv < 0.05
cv_filtered <- cv[sig_idx]
rate_sums_filtered <- rate_sums[sig_idx]

hist(cv_filtered)

plot(rate_sums_filtered, abs(cv_filtered))
cor.test(rate_sums_filtered, abs(cv_filtered), method = "pearson")
cor.test(rate_sums_filtered, abs(cv_filtered), method = "spearman")


#----------  Graph

source('adj-matrix-ST.R')
adj_matrix

library(igraph)
# Create undirected graph from adjacency
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
plot(g)
# Get shortest path distance matrix
dist_matrix <- distances(g)


pv <- pval(corr.pearson.log, use.bonf = T) 
cv <- corr_val(corr.pearson.log)

# sig_idx <- pv < 0.05
# cv_filtered <- cv[sig_idx]
# 
# head(cv)
# head(adj_matrix)

# Load required package
library(vegan)

# Inputs
# - adj_matrix: square adjacency matrix with named rows/columns (body parts)
# - cv: named numeric vector of pairwise correlations (e.g., "cranium vs mouthparts")
# - pv: named numeric vector of p-values matching cv

# ---- Step 1: Filter correlations by p-value
sig_idx <- pv < 0.05
cv_filtered <- cv[sig_idx]

# ---- Step 2: Convert adjacency matrix to distance matrix
# Replace 1s with 1 and 0s with a large number to simulate disconnected pairs
adj_dist <- as.matrix(adj_matrix)
adj_dist[adj_dist == 1] <- 1
adj_dist[adj_dist == 0] <- 100  # large pseudo-distance for non-adjacent parts
diag(adj_dist) <- 0
adj_dist <- as.dist(adj_dist)

# ---- Step 3: Create correlation distance matrix (1 - abs(correlation))
# Extract the part names from the names of cv_filtered
extract_parts <- function(name) strsplit(name, " vs ")[[1]]
part_pairs <- lapply(names(cv), extract_parts)
parts <- unique(unlist(part_pairs))

# Initialize correlation distance matrix
cor_dist_mat <- matrix(0, nrow = length(parts), ncol = length(parts),
                       dimnames = list(parts, parts))

# Fill the matrix with 1 - abs(correlation)
for (i in seq_along(cv_filtered)) {
  pair <- part_pairs[[i]]
  val <- 1 - abs(cv_filtered[i])
  cor_dist_mat[pair[1], pair[2]] <- val
  cor_dist_mat[pair[2], pair[1]] <- val
}
diag(cor_dist_mat) <- 0
dist_matrix_dist

# Convert to 'dist' object
cor_dist <- as.dist(cor_dist_mat)
dist_matrix_dist <- as.dist(dist_matrix)
# ---- Step 4: Ensure common names between the two matrices
# common_names <- intersect(attr(adj_dist, "Labels"), attr(cor_dist, "Labels"))
# adj_dist <- as.dist(as.matrix(adj_dist)[common_names, common_names])
# cor_dist <- as.dist(as.matrix(cor_dist)[common_names, common_names])


# ---- Step 5: Run Mantel test
plot(cor_dist, dist_matrix_dist)

mantel_result <- mantel(cor_dist, dist_matrix_dist, method = "spearman", permutations = 999)

# ---- Step 6: Print result
print(mantel_result)


