library(corHMM)
library(phytools)
source('R/utils-maddfitz-root.r')

# ------------------ Simulate a random tree with 50 tips
tree <- pbtree(n = 50, scale = 1)
plot(tree)
Q <- matrix(c(-1, 1,
              1, -1), 2, 2)
char1 <- sim.history(tree, Q)$states
# ------------------

# ------------------ FIT
taxa <- cbind(names(char1), char1)
# for experimental purposes let's make one ambiguous
taxa[1,2] <- "1&2"
recon <- rayDISC(tree, taxa, model='ARD', node.states="none", root.p='maddfitz')

# ------------------ Get Root
# This is from corHMM to get params
#  one way to get the parameters from your corHMM object in the correct order
p <- sapply(1:max(recon$index.mat, na.rm = TRUE), function(x) 
  na.omit(c(recon$solution))[na.omit(c(recon$index.mat) == x)][1])

# REMEBER to use the same model e.g. model='ARD'
# Use method = "scaled"
anc <- ancRECON_ROOT(phy = tree, data = taxa, p=p, model='ARD', method = "scaled", 
         rate.cat = 1, root.p = 'maddfitz', get.likelihood=T, get.tip.states = T)

# Maddfitz root
anc$root

# tip state prior
anc$lik.tip.states
# tips state prior should be scaled for phytools (I believe)
prior <- anc$lik.tip.states
prior <- apply(prior, 1, function(x) x/sum(x))
prior <- t(prior)
prior
