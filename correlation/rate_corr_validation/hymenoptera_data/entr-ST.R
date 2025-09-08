# Load libraries.
library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)
source("R/utils-ST.R")

# # Import some functions.
# source("R/helpers.R")
# 
# # Load Hymenoptera dataset.
# load("data_hym/paramo_stm_final.RDA")


#stm_merg <- lapply(stm_amalg, function(x) merge_tree_cat_list(x) )
#saveRDS(stm_merg, "data_hym/stm_merg.RDS")
stm_merg <-readRDS("data_hym/stm_merg.RDS")
tree <-readRDS("data_hym/hym_tree.RDS")
plot(tree)

# branch rate
rates_list <- lapply(stm_merg, get_branch_rate_across_smaps)
rates <- do.call(cbind, rates_list)
head(rates)
rate_total <- apply(rates, 2, sum)
rate_total

rate_total_branch <- apply(rates, 1, sum)


# Min-max scale each column to [0, 1]
rates_scaled <- apply(rates, 2, function(col) {
  (col - min(col)) / (max(col) - min(col))
})
head(rates_scaled)
max(rates_scaled[,1])


ncol(rates_scaled)
entropy(rep(2,15))
entropy(c(1,rep(0.01,14)))

# Apply entropy row-wise
row_entropies <- apply(rates_scaled, 1, entropy)
hist(row_entropies)

# using absolute rates
row_entropies <- apply(rates, 1, entropy)
hist(row_entropies)
plot(row_entropies, rate_total_branch)
plot(log(row_entropies), log(rate_total_branch))
cor.test(row_entropies, rate_total_branch, method = "spearman")
cor.test(row_entropies, rate_total_branch, method = "pearson")

library(viridis)
#hist(rate_total_branch)
#row_entropies <- tree$edge.length
#row_entropies <- log(rate_total_branch)
# Prepare colors


plot_branch_colored_tree(tree, row_entropies, palette = viridis(4), 
                         title = "Over rep", legend_title = 'Cols')

#--------- Enrichment
hist(rates_scaled, breaks=30)
rates_scaled
rate_total
br1 <- rates_scaled[1,]
hist(rates_scaled[,1])
head(rates_scaled)
head(rates)

enriched_over <- get_enrichment_mask(rates, lower_q = 0.1, upper_q = 0.9, direction = "over", return_type = "binary")
br_over <- apply(enriched_over, 1, sum)
apply(enriched_over, 2, sum)
head(enriched_over)

br_over_binary <- br_over
br_over_binary <-(br_over_binary>0)*1

plot(row_entropies, br_over)

plot_branch_colored_tree(tree, br_over , palette = viridis(15), 
                         title = "Over rep", legend_title = 'Cols')

plot_branch_colored_tree(tree, br_over_binary , palette = viridis(3), 
                         title = "Over rep", legend_title = 'Cols')

plot_branch_colored_tree(tree, enriched_over[,12], palette = viridis(100), 
                         title = "Over rep", legend_title = 'Cols')

# head
plot_branch_colored_tree(tree, (apply(enriched_over[,c(1:2)], 1, sum) > 0)*1, palette = viridis(3), 
                         title = "Head", legend_title = 'Cols')

# thorax
plot_branch_colored_tree(tree, (apply(enriched_over[,c(3:8)], 1, sum) > 0)*1, palette = viridis(3), 
                         title = "Thorax", legend_title = 'Cols')

# Meta+Genit
plot_branch_colored_tree(tree, (apply(enriched_over[,c(9:10)], 1, sum) > 0)*1, palette = viridis(3), 
                         title = "Thorax", legend_title = 'Cols')

#others
plot_branch_colored_tree(tree, (apply(enriched_over[,c(11:15)], 1, sum) > 0)*1, palette = viridis(3), 
                         title = "Others", legend_title = 'Cols')

#----
enriched_under <- get_enrichment_mask(rates, direction = "under", return_type = "binary")
br_under <- apply(enriched_under, 1, sum)

br_under_binary <- br_under
br_under_binary <-(br_under_binary>0)*1

plot(row_entropies, br_under)

plot_branch_colored_tree(tree, br_under , palette = viridis(100), 
                         title = "Unde rep", legend_title = 'Cols')

plot_branch_colored_tree(tree, br_under_binary , palette = viridis(100), 
                         title = "Under rep", legend_title = 'Cols')


#----------- Clustering
head(enriched_over)
enriched_over_masked <- get_enrichment_mask(rates, lower_q = 0.1, upper_q = 0.9, direction = "over", return_type = "masked")
# Transpose the matrix so traits are rows
# trait_matrix <- t(enriched_over)
trait_matrix <- enriched_over
#trait_matrix <- t(enriched_over_masked)

# Compute binary distance (Jaccard by default)
dist_traits <- dist(trait_matrix, method = "binary")

# Hierarchical clustering
hc_traits <- hclust(dist_traits, method = "average")

# Plot dendrogram
plot(hc_traits, main = "Clustering of Enriched Body Parts", xlab = "", sub = "")


library(igraph)

# Compute co-enrichment matrix
co_occurrence <- trait_matrix %*% t(trait_matrix)

# Remove self-links if needed
diag(co_occurrence) <- 0

# Create graph from adjacency matrix
g <- graph_from_adjacency_matrix(co_occurrence, mode = "undirected", weighted = TRUE)

# Optional: filter low-weight edges
#g <- delete_edges(g, E(g)[weight < 2])  # Only keep edges with 2+ co-enrichments

# Plot network
plot(g,
     edge.width = E(g)$weight,
     vertex.label.cex = 0.9,
     main = "Trait Co-Enrichment Network")

# Community detection using Louvain method
communities <- cluster_louvain(g)

# Print communities
print(communities)

# Get membership for each vertex
membership <- membership(communities)

# Assign colors to vertices based on community
V(g)$color <- membership

# Plot with community colors
plot(g,
     layout = layout_with_fr(g),
     vertex.color = V(g)$color,
     vertex.label.cex = 0.9,
     edge.width = E(g)$weight,
     main = "Trait Co-Enrichment Network with Communities")


#-------- PCA

enriched_over_nonempty <- enriched_over[rowSums(enriched_over) > 0, ]
# Compute Jaccard distance
library(vegan)
jaccard_dist <- vegdist(enriched_over_nonempty, method = "jaccard")

# Hierarchical clustering
hc <- hclust(jaccard_dist)
plot(hc, main = "Hierarchical Clustering of Rows (Jaccard)")

# PCA

pca_res <- prcomp(enriched_over_masked, scale. = TRUE)
summary(pca_res)

# Basic PCA scatter plot (first two PCs)
plot(pca_res$x[,1], pca_res$x[,2],
     xlab = "PC1", ylab = "PC2",
     main = "PCA of Enriched Over Matrix",
     pch = 19, col = "steelblue")

text(pca_res$x[,1], pca_res$x[,2],
     labels = 1:nrow(enriched_over_masked),
     pos = 3, cex = 0.7)
