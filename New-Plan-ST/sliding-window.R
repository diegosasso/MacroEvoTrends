library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)

library(igraph)
library(vegan)
source("R/utils-ST.R")
source('R/region_class.R')

# ====================== Read Data ====================== #
# ---- 15 body regions (first order) ----
stm_merg <- readRDS("data/stm_merg.RDS")
tree <-readRDS("data/hym_tree.RDS")
# branch rate
rates_list <- lapply(stm_merg, get_branch_rate_across_smaps)
rates <- do.call(cbind, rates_list)
head(rates)
# proportion of zeros
colMeans(rates == 0)
# total rate per body regions (is it correct ????)
rate_total <- apply(rates, 2, sum) #/ nrow(rates) # should it be devided by N branches

# ---- Number of states per BR ----
# load('data/char_num.RDA')
# #anat_ent_state
# # log of states per character per BR
# log_states <- lapply(anat_ent_state, function(x) sum(log(get_len(x))) ) %>% unlist
# log_states
# tot_states <- lapply(anat_ent_state, function(x) sum(get_len(x)) ) %>% unlist
# tot_states
# 
# rate_total_norm <- rate_total/tot_states
# rate_total_norm

tree
head(rates)


library(dplyr)
library(tidyr)

library(ggplot2)

# edge_times from earlier function
edge_times <- edge_times_from_tips(tree)
head(rates)
edge_times


count_branches_sliding(edge_times, window_size = 40, step_size = 1, plot = TRUE)

# Run rolling window
cor_results <- rolling_rate_correlation(rates[, c(1:15)], edge_times,
                                        window_size = 40, 
                                        step_size = 1,
                                        trim_quantiles = NULL, # c(0.01, 0.99), 
                                        log_transform = TRUE,
                                        use.abs.value=TRUE)

head(cor_results)
cor_list_to_df <- function(cor_results) {
  df <- do.call(rbind, lapply(names(cor_results), function(t) {
    mat <- cor_results[[t]]
    if (all(is.na(mat))) return(NULL)
    as.data.frame(as.table(mat)) %>%
      rename(region1 = Var1, region2 = Var2, cor = Freq) %>%
      mutate(time = as.numeric(t))
  }))
  return(df)
}

df_cor <- cor_list_to_df(cor_results)
head(df_cor)
tail(df_cor)

xx <- df_cor %>% filter(time==10)
hist(xx$cor)

# Average correlation over all pairs
# df_avg <- df_cor %>%
#   group_by(time) %>%
#   summarize(mean_cor = mean(cor, na.rm = TRUE),
#             se = sd(cor, na.rm = TRUE)/sqrt(n()))

# ggplot(df_avg, aes(x = time, y = mean_cor)) +
#   geom_line(color = "black") +
#   geom_ribbon(aes(ymin = mean_cor - se, ymax = mean_cor + se),
#               alpha = 0.2, fill = "grey") +
#   scale_x_reverse() +
#   theme_classic() +
#   labs(y = "Mean pairwise correlation", x = "Time (Ma)")

# ====================== test  ====================== #


# ======================  ====================== #
df_avg <- df_cor %>%
  group_by(time) %>%
  summarize(
    mean_cor = mean(cor, na.rm = TRUE),
    se = sd(cor, na.rm = TRUE) / sqrt(n()),
    sd = sd(cor, na.rm = TRUE)
  )

ggplot(df_avg, aes(x = time, y = mean_cor)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = mean_cor - se, ymax = mean_cor + se),
              alpha = 0.2, fill = "grey") +
  scale_x_reverse(
    breaks = seq(0, 300, by = 25),   # ticks every 25 units
    expand = c(0, 0)
  ) +
  theme_classic() +
  labs(y = "Mean pairwise correlation", x = "Time (Ma)")


ggplot(df_cor, aes(x = factor(time), y = cor)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey80", color = "black") +
  scale_x_discrete(labels = function(x) paste0(x, " Ma")) +
  coord_flip() + # optional: flip for readability
  theme_classic(base_size = 14) +
  labs(y = "Pairwise correlation", x = "Time slice")

#------------------------------
trend_results <- df_cor %>%
  group_by(region1, region2) %>%
  do({
    fit <- lm(cor ~ time, data = .)
    data.frame(
      slope = coef(fit)[["time"]],
      p_value = summary(fit)$coefficients["time", "Pr(>|t|)"]
    )
  })

trend_results

library(mgcv)
gam_fit <- gam(mean_cor ~ s(time), data = df_avg)
summary(gam_fit)



# Fit GAM
gam_fit <- gam(mean_cor ~ s(time), data = df_avg)

# Predict from GAM with CI
pred <- predict(gam_fit, newdata = data.frame(time = df_avg$time),
                se.fit = TRUE)

df_pred <- data.frame(
  time = df_avg$time,
  fit = pred$fit,
  lower = pred$fit - 2 * pred$se.fit,  # ~95% CI
  upper = pred$fit + 2 * pred$se.fit
)

# Plot
ggplot() +
  # raw rolling correlations
  geom_point(data = df_avg, aes(x = time, y = mean_cor),
             color = "grey50", size = 2, alpha = 0.6) +
  
  # GAM confidence ribbon
  geom_ribbon(data = df_pred, aes(x = time, ymin = lower, ymax = upper),
              fill = "#a6cee3", alpha = 0.4) +
  
  # GAM smooth line
  geom_line(data = df_pred, aes(x = time, y = fit),
            color = "#1f78b4", size = 1.2) +
  
  scale_x_reverse(expand = c(0,0)) +   # time decreasing to the right
  labs(x = "Time from root (Ma)",
       y = "Mean pairwise correlation") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  )




#------ Major rewiring

early <- df_cor %>% filter(time > 225)   # example cutoff
late  <- df_cor %>% filter(time > 100 & time < 150)

cor_early <- early %>%
  group_by(region1, region2) %>%
  summarize(cor = mean(cor, na.rm = TRUE))

cor_late <- late %>%
  group_by(region1, region2) %>%
  summarize(cor = mean(cor, na.rm = TRUE))

#------ Mantel
# Make correlation matrices
mat_early <- cor_early %>%
  pivot_wider(names_from = region2, values_from = cor) %>%
  tibble::column_to_rownames("region1") %>%
  as.matrix()

mat_late <- cor_late %>%
  pivot_wider(names_from = region2, values_from = cor) %>%
  tibble::column_to_rownames("region1") %>%
  as.matrix()

# Ensure both have same row/col order
common_regions <- intersect(rownames(mat_early), rownames(mat_late))
mat_early <- mat_early[common_regions, common_regions]
mat_late  <- mat_late[common_regions, common_regions]

# Convert to distance objects (e.g., 1 - correlation as dissimilarity)
dist_early <- as.dist(1 - mat_early)
dist_late  <- as.dist(1 - mat_late)

plot(mat_early, mat_late)
# Mantel test
library(vegan)
plot(dist_early, dist_late)
mantel(dist_early, dist_late, method = "pearson", permutations = 999)


cor_change <- merge(cor_early, cor_late,
                    by = c("region1", "region2"),
                    suffixes = c("_early", "_late"))

cor_change <- cor_change %>%
  mutate(diff = cor_late - cor_early) %>%
  arrange(desc(abs(diff)))
cor_change


#------ Per Body region correlation

cor_brs <- cor_trimmed(
  rates,
  method = "pearson",
  log_transform = TRUE,
  log_const = 1e-6,
  trim_quantiles = NULL, #c(0.01, 0.99),
  use = "pairwise.complete.obs"
)
hist(cor_brs)

#cor_mat <- cor_brs
test_within_vs_across <- function(cor_mat, region_class) {
  results <- list()
  
  # Ensure row/col order matches region_class
  region_class <- region_class[match(colnames(cor_mat), region_class$region), ]
  
  #g='head'
  for (g in unique(region_class$group)) {
    # Regions in this group
    in_group <- region_class$region[region_class$group == g]
    idx_in <- which(colnames(cor_mat) %in% in_group)
    idx_out <- setdiff(seq_len(ncol(cor_mat)), idx_in)
    
    # Extract within-group correlations
    within_vals <- cor_mat[idx_in, idx_in][upper.tri(cor_mat[idx_in, idx_in])]
    
    # Extract across-group correlations
    across_vals <- cor_mat[idx_in, idx_out]
    
    # Store summary
    results[[g]] <- data.frame(
      group = g,
      mean_within = mean(within_vals, na.rm = TRUE),
      mean_across = mean(across_vals, na.rm = TRUE),
      n_within = length(within_vals),
      n_across = length(across_vals),
      t_test_p = tryCatch(
        t.test(within_vals, across_vals)$p.value,
        error = function(e) NA
      )
    )
  }
  
  # Bind into one data frame
  do.call(rbind, results)
}

# Example usage:
out <- test_within_vs_across(cor_brs, region_class)
print(out)


library(ggplot2)
library(reshape2)

cor_brs <- abs(cor_brs)
cor_brs[cor_brs<0.4] <- 0

# Convert correlation matrix into long format
cor_long <- melt(cor_brs)

ggplot(cor_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "white", mid = "#2166ac", high = "#b2182b",
                       midpoint = 0.5, limit = c(0, 1), 
                       name = "Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "")



# ====================== Modularity ====================== #
library(igraph)

# Function to compute modularity over time
modularity_over_time <- function(cor_results) {
  mod_list <- lapply(names(cor_results), function(t) {
    mat <- cor_results[[t]]
    
    # Convert to igraph, keep positive weights only
    g <- graph_from_adjacency_matrix(
      mat, mode = "undirected", weighted = TRUE, diag = FALSE
    )
    
    # Run Louvain clustering
    comm <- cluster_louvain(g, weights = E(g)$weight)
    
    data.frame(
      time = as.numeric(t),
      modularity = modularity(comm),
      n_clusters=length(comm)
    )
  })
  
  df_mod <- do.call(rbind, mod_list)
  
  # Plot with time reversed
  p <- ggplot(df_mod, aes(x = time, y = modularity)) +
    geom_line(size = 1, color = "steelblue") +
    geom_point(color = "steelblue") +
    scale_x_reverse(expand = c(0, 0)) +
    theme_classic(base_size = 14) +
    labs(x = "Time (Ma)", y = "Modularity (Q)")
  
  return(list(data = df_mod, plot = p))
}

# Example usage
res <- modularity_over_time(cor_results)
print(res$plot)
res$data

plot(res$data$time, res$data$n_clusters, type='l')

plot(df_avg$mean_cor, res$data$modularity[1:14])
abline(lm(res$data$modularity[1:nrow(df_avg)] ~ df_avg$mean_cor),
       col = "red", lwd = 2)
summary(lm(res$data$modularity[1:nrow(df_avg)] ~ df_avg$mean_cor))


# merge the two datasets by time
df_plot <- merge(df_avg, res$data, by = "time", all = TRUE)

ggplot(df_plot, aes(x = time)) +
  # Mean correlation with SE ribbon
  geom_ribbon(aes(ymin = mean_cor - se, ymax = mean_cor + se),
              alpha = 0.2, fill = "grey") +
  geom_line(aes(y = mean_cor), color = "black", size = 1) +
  
  # Modularity line, mapped to secondary axis
  geom_line(aes(y = modularity * (max(mean_cor, na.rm = TRUE) / 
                                    max(modularity, na.rm = TRUE))),
            color = "firebrick", size = 1, linetype = "solid") +
  
  scale_x_reverse(
    breaks = seq(0, 300, by = 25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Mean pairwise correlation",
    sec.axis = sec_axis(~ . * (max(res$data$modularity, na.rm = TRUE) / 
                                 max(df_avg$mean_cor, na.rm = TRUE)),
                        name = "Modularity (Q)")
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "firebrick")
  ) +
  labs(x = "Time (Ma)")

# ====================== Modularity OLD ====================== #

# Build graph from correlation matrix
#g <- graph_from_adjacency_matrix(cor_brs, mode = "undirected", weighted = TRUE, diag = FALSE)

g <- graph_from_adjacency_matrix(cor_results[[100]], mode = "undirected", weighted = TRUE, diag = FALSE)

# Optional: remove weak correlations
E(g)$weight[abs(E(g)$weight) < 0.2] <- 0
g <- delete_edges(g, E(g)[weight == 0])

# Louvain community detection
comm <- cluster_louvain(g, weights = E(g)$weight)
comm <- cluster_walktrap(g, weights = E(g)$weight)
comm <- cluster_fast_greedy(g, weights = E(g)$weight)
comm <- cluster_leading_eigen(g, weights = E(g)$weight)

# Module membership
membership(comm)


plot(comm, g, vertex.size = 15, 
     vertex.label.cex = 0.8,
     edge.width = abs(E(g)$weight)*5,
     edge.color = ifelse(E(g)$weight > 0, "firebrick", "steelblue"))

#---------------- Compare two clusterings
g <- graph_from_adjacency_matrix(mat_early, mode = "undirected", weighted = TRUE, diag = FALSE)
E(g)$weight[abs(E(g)$weight) < 0.4] <- 0
g <- delete_edges(g, E(g)[weight == 0])
comm <- cluster_leading_eigen(g, weights = E(g)$weight)
modularity(comm)
membership_vec <- membership(comm)
modularity(g, membership_vec, weights = E(g)$weight)
plot(comm, g, vertex.size = 15, 
     vertex.label.cex = 0.8,
     edge.width = abs(E(g)$weight)*5,
     edge.color = ifelse(E(g)$weight > 0, "firebrick", "steelblue"))

g <- graph_from_adjacency_matrix(mat_late, mode = "undirected", weighted = TRUE, diag = FALSE)
E(g)$weight[abs(E(g)$weight) < 0.4] <- 0
g <- delete_edges(g, E(g)[weight == 0])
comm1 <- cluster_leading_eigen(g, weights = E(g)$weight)
modularity(comm1)
modularity(g, membership_vec, weights = E(g)$weight)
plot(comm1, g, vertex.size = 15, 
     vertex.label.cex = 0.8,
     edge.width = abs(E(g)$weight)*5,
     edge.color = ifelse(E(g)$weight > 0, "firebrick", "steelblue"))

compare(comm, comm1, method = "nmi")    # normalized mutual information
compare(comm, comm1, method = "vi")     # variation of information
compare(comm, comm1, method = "rand")   # Rand index



m1 <- membership(comm)   # first partition
m2 <- membership(comm1)  # second partition

tab <- table(m1, m2)
print(tab)
chisq.test(tab)
chisq.test(tab, simulate.p.value = TRUE, B = 10000)

library(ggalluvial)
df <- data.frame(
  comm  = paste0("comm:", membership(comm)),
  comm1 = paste0("comm1:", membership(comm1))
)
ggplot(df, aes(axis1 = comm, axis2 = comm1)) +
  geom_alluvium(aes(fill = comm)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  labs(x = "", y = "Vertices")



# install.packages("mcclust") # if not installed
library(mcclust)
variation_of_information <- function(cl1, cl2) {
  tab <- table(cl1, cl2)
  n <- sum(tab)
  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)
  
  Hx <- -sum(px * log(px + 1e-12))
  Hy <- -sum(py * log(py + 1e-12))
  Ixy <- sum(pxy * log((pxy + 1e-12) / (outer(px, py) + 1e-12)))
  
  VI <- Hx + Hy - 2 * Ixy
  return(VI)
}

# m1 and m2 are cluster memberships (vectors of same length)
obs_vi <- variation_of_information(m1, m2)
print(obs_vi)

perm_vi <- replicate(1000, {
  m2_perm <- sample(m2)  # permute cluster assignments
  names(m2_perm) <- NULL
  variation_of_information(m1, m2_perm)
})

p_value <- mean(perm_vi <= obs_vi)  # smaller VI means more similar
p_value
