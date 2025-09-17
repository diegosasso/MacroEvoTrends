library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)
source("R/utils-ST.R")
source("R/sim-test.R")

# False positive rate (expect ≈ 0.05)
null_results <- run_simulation(n.reps = 100, true_corr = 0, cor_method='pearson')
cat("False positive rate:", null_results$proportion_significant, "\n")

# Power test (expect higher value)
power_results <- run_simulation(n.reps = 100, true_corr = 0.2, , cor_method='pearson')
cat("Power (true corr = 0.5):", power_results$proportion_significant, "\n")


# False positive rate (expect ≈ 0.05)
null_results <- run_simulation(n.reps = 100, true_corr = 0, cor_method='spearman')
cat("False positive rate:", null_results$proportion_significant, "\n")

# Power test (expect higher value)
power_results <- run_simulation(n.reps = 100, true_corr = 0.2, , cor_method='spearman')
cat("Power (true corr = 0.5):", power_results$proportion_significant, "\n")


#---------------

n_reps <- 1000      # Number of replicates
n_points <- 150     # Sample size
alpha <- 0.05       # Significance level

# ---- 1. Simulate under null hypothesis (no correlation) ----
false_positives <- replicate(n_reps, {
  x <- rnorm(n_points)
  y <- rnorm(n_points)  # Independent from x
  pval <- cor.test(x, y, method = "pearson")$p.value
  pval < alpha
})

fpr <- mean(false_positives)  # False positive rate
cat("False Positive Rate (null case):", fpr, "\n")

# ---- 2. Simulate under alternative hypothesis (true correlation) ----
power_hits <- replicate(n_reps, {
  x <- rnorm(n_points)
  y <- 0.2 * x + sqrt(1 - 0.2^2) * rnorm(n_points)  # y correlated with x, rho ≈ 0.5
  pval <- cor.test(x, y, method = "pearson")$p.value
  pval < alpha
})

power <- mean(power_hits)  # Power of the test
cat("Power (rho = 0.5):", power, "\n")

