# Load libraries
library(phytools)
source('R/utils-branch-ama.R')

# 1. Simulate a random tree with 50 tips
tree <- pbtree(n = 10, scale = 1)
plot(tree)

# 2. Define a transition rate matrix (Q) for a binary trait
Q <- matrix(c(-1, 1,
              1, -1), 2, 2)

# 3. Simulate three discrete characters
char1 <- sim.history(tree, Q)$states
char2 <- sim.history(tree, Q)$states
char3 <- sim.history(tree, Q)$states


# 4. Fit a model for each character and generate stochastic maps
maps1 <- make.simmap(tree, char1, model = "SYM", nsim = 5)
maps2 <- make.simmap(tree, char2, model = "SYM", nsim = 5)
maps3 <- make.simmap(tree, char3, model = "SYM", nsim = 5)

plot(maps1[[2]])
edgelabels()

# Cranium with 3 chars each char with 5 maps
cranium <- list(char1=maps1, char2=maps2, char3=maps3)
cranium

# This gives a sample of amalgamated rates across characters
cranium.branch.sample <- ama_bracnhes_multi_char(cranium)
cranium.branch.sample

# This amalgamates and summarizes rates across characters
mean.cranium <- summarize_branch_rates_multi(cranium, 'mean')

# This summarizes rates for one individual character
mu1 <- summarize_branch_rates(cranium$char1, 'mean')
mu2 <- summarize_branch_rates(cranium$char2, 'mean')
mu3 <- summarize_branch_rates(cranium$char3, 'mean')

mu1+mu2+mu3
mean.cranium



