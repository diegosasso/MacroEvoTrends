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
#    We'll use make.simmap() with 100 replicates
maps1 <- make.simmap(tree, char1, model = "SYM", nsim = 50)
maps2 <- make.simmap(tree, char2, model = "SYM", nsim = 50)
maps3 <- make.simmap(tree, char3, model = "SYM", nsim = 50)

plot(maps1[[2]])
edgelabels()

cranium <- list(char1=maps1, char2=maps2, char3=maps3)
cranium

#list.of.chars <- cranium
cranium.branch.sample <- ama_bracnhes_multi_char(cranium)
cranium.branch.sample
summarize_branch_rates_multi(cranium, 'mean')

m1 <- summarize_branch_rates(cranium$char1, 'mean')
m2 <- summarize_branch_rates(cranium$char2, 'mean')
m3 <- summarize_branch_rates(cranium$char3, 'mean')

m1+m2+m3
summarize_branch_rates_multi(cranium, 'mean')



m1 <- summarize_branch_rates(cranium$char1, 'median')
m2 <- summarize_branch_rates(cranium$char2, 'median')
m3 <- summarize_branch_rates(cranium$char3, 'median')

m1+m2+m3
summarize_branch_rates_multi(cranium, 'median')
