library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)

library(igraph)
library(vegan)
source("R/utils-ST.R")


# ====================== Read Data ====================== #
# ---- 15 body regions (first order) ----
stm_merg <- readRDS("data/stm_merg.RDS")
# branch rate
rates_list <- lapply(stm_merg, get_branch_rate_across_smaps)
rates <- do.call(cbind, rates_list)
head(rates)
# proportion of zeros
colMeans(rates == 0)
# total rate per body regions (is it correct ????)
rate_total <- apply(rates, 2, sum) / nrow(rates) # should it be devided by N branches

# ---- Number of states per BR ----
load('data/char_num.RDA')
#anat_ent_state
# log of states per character per BR
log_states <- lapply(anat_ent_state, function(x) sum(log(get_len(x))) ) %>% unlist
log_states
tot_states <- lapply(anat_ent_state, function(x) sum(get_len(x)) ) %>% unlist
tot_states

rate_total_norm <- rate_total/tot_states
rate_total_norm

# ====================== Pairwise correlation: 15 body regions (first order) ====================== #
corr.pearson.log <-run_all_pairwise_correlations_parallel(
  stm_merg,
  cor_method = "pearson",
  remove_zeros = FALSE,
  scale = FALSE,
  log.trans = T,
  log.constant = 1e-6,
  quantile.trim = c(0.01, 0.99),
  use.resampling = F,
  n.cores = parallel::detectCores() - 1
)



# ====================== Mantel test ====================== #
source('adj-matrix.R')

# ---- H1 ----
adj_matrix
# Create undirected graph from adjacency
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
plot(g)
dist_matrix <- distances(g)
dist_matrix_dist <- as.dist(dist_matrix)


pv <- pval(corr.pearson.log, use.bonf = T) 
cv <- corr_val(corr.pearson.log)
pv[] <- .001
# get inverse cv to represent it as a dissimilarity matrix
cv.inv <- 1-abs(cv)
M <- create_pairwise_matrix(cv.inv, pv, dist_matrix_dist, mode = "na")
head(M)
#create_pairwise_matrix(cv.inv, pv, dist_matrix_dist, mode = 'none')

plot(M$cor_mat_dist, M$dist_mat_dist)
mantel_h1 <- vegan::mantel(M$cor_mat_dist, M$dist_mat_dist, method = "pearson", permutations = 999, na.rm = T)
mantel_h1



# ---- H2 ----
plot(g2)
dist_matrix <- distances(g2)
d.h2 <- dist_matrix
dist_matrix_dist <- as.dist(dist_matrix)


pv <- pval(corr.pearson.log, use.bonf = T) 
pv[] <- .001
cv <- corr_val(corr.pearson.log)

# get inverse cv for dissimilarity matrix
cv.inv <- 1-abs(cv)
M <- create_pairwise_matrix(cv.inv, pv, dist_matrix_dist, mode = "na")
head(M)

plot(M$cor_mat_dist, M$dist_mat_dist)

mantel_h2 <- vegan::mantel(M$cor_mat_dist, M$dist_mat_dist, method = "pearson", permutations = 999, na.rm = T)
mantel_h2

# ---- Summary ----
mantel_h1
mantel_h2

# ====================== Regression ====================== #

# ---- Make. data ----
pv <- pval(corr.pearson.log, use.bonf = T) 
cv <- corr_val(corr.pearson.log)
ana.d
plot(g)
dist_matrix <- distances(g)

h1.pairs <- extract_pairs(dist_matrix, names(ana.d))
cv_abs <- 1-abs(cv[names(ana.d)])


df <- data.frame(
  cv_abs = cv_abs,
  h1 = h1.pairs,
  h2 = ana.d
)
df

# ---- Fir models ----

# Linear models
lm_h1 <- lm(cv_abs ~ h1, data = df)
lm_h2 <- lm(cv_abs ~ h2, data = df)
# 2. Multiple regression
lm_both <- lm(cv_abs ~ h1 + h2, data = df)

# ---- Summary ----

summary(lm_h1)
summary(lm_h2)
summary(lm_both)

AIC(lm_h1, lm_h2, lm_both)



#---- Test multicollinearity with the Variance Inflation Factor (VIF).
library(car)
# VIF for the multiple regression
vif(lm_both)
# are very close to 1, which means there is almost no multicollinearity between h1 and h2.




# ====================== Regression: Rates. vs Body Categories ====================== #



# ---- Plot Rates ----
barplot(rate_total,
        las = 2,                      # Rotate labels vertically
        col = "skyblue",
        main = "Total Rates by Body Part",
        ylab = "Rate")

barplot(log(rate_total),
        las = 2,                      # Rotate labels vertically
        col = "skyblue",
        main = "Total Rates by Body Part",
        ylab = "Rate")

barplot(rate_total_norm,
        las = 2,                      # Rotate labels vertically
        col = "skyblue",
        main = "Total Rates by Body Part",
        ylab = "Rate")

identical(names(log_states), names(rate_total))
rate_log_norm <- log(rate_total)-log_states
rate_norm <- log(rate_total)-log(tot_states)

plot(rate_total, tot_states, ylim=c(0,150))
plot(log(rate_total), log(tot_states), ylim=c(0,150))
plot(log(rate_total), log_states, ylim=c(0,150))
#abline(0, 1, col = "red", lwd = 2, lty = 2)

barplot(rate_log_norm,
        las = 2,                      # Rotate labels vertically
        col = "skyblue",
        main = "Total Rates by Body Part",
        ylab = "Rate")

barplot(exp(rate_norm),
        las = 2,                      # Rotate labels vertically
        col = "skyblue",
        main = "Total Rates by Body Part",
        ylab = "Rate")

bd <- read.csv('data/dist-Hym - Sheet2.csv')
#df_rates <- data.frame(region = names(rate_total), rate_total = rate_total)
#df_rates <- data.frame(region = names(rate_total), rate_total = log(rate_total))
df_rates <- data.frame(region = names(rate_total), rate_total = rate_total_norm)
#df_rates <- data.frame(region = names(rate_total), rate_total = rate_norm)
df_rates <- merge(df_rates, bd, by = "region")
df_rates

# Convert grouping columns h1-h6 to factors
df_rates[, c("h1", "h2", "h3", "h5", "h6", "h7", "h8")] <- 
  lapply(df_rates[, c("h1", "h2", "h3", "h5", "h6", "h7", "h8")], factor)
df_rates

# Fit separate models for each hypothesis
lm_h1 <- lm(rate_total ~ factor(h1), data = df_rates)
lm_h2 <- lm(rate_total ~ factor(h2), data = df_rates)
lm_h3 <- lm(rate_total ~ factor(h3), data = df_rates)
lm_h5 <- lm(rate_total ~ factor(h5), data = df_rates)
lm_h6 <- lm(rate_total ~ factor(h6), data = df_rates)
#lm_h7 <- lm(rate_total ~ factor(h7), data = df_rates)
lm_h8 <- lm(rate_total ~ factor(h8), data = df_rates)


# Compare models with AIC
AIC(lm_h1, lm_h2, lm_h3, lm_h5, lm_h6, lm_h8)

# Summaries (to check within-group means and effect size)
summary(lm_h1)
summary(lm_h2)
summary(lm_h3)
summary(lm_h6)
summary(lm_h2)
summary(lm_h8)

par(mfrow = c(2,2))
plot(lm_h3)  
# Null hypothesis: residuals are normal.	p > 0.05 → normality not rejected.
shapiro.test(residuals(lm_h3))
# Null hypothesis: variance is constant. p > 0.05 → homoscedasticity not rejected.
lmtest::bptest(lm_h3)


# Group effect
anova(lm_h8)  
