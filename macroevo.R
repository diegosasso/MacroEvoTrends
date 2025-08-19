
# Load packages.
library(phytools)
library(rphenoscate)
library(phytools)

##--------------------------------------------------------------------------------##

# Load data.
hym_tree <- readRDS("data/hym_tree.RDS")

# Paint a clade (Proctotrupomorpha) with given macroevolutionary regime. 
hym_tree_reg <- paintSubTree(hym_tree, node = 86, state = "Q2", anc.state = "Q1",stem = T)

# Check regimes.
plotSimmap(hym_tree_reg, lwd = 3, fsize = 0.5)
dev.off()

# Extract subtree.
hym_sub <- extract.clade(hym_tree, node = 86)

# Add root edge.
hym_sub$root.edge <- hym_tree$edge.length[8]
hym_sub <- rootedge.to.singleton(hym_sub)

# Add temporary outgroup (just for stochastic mapping).
hym_sub <- bind.tip(hym_sub, where = 35, "outgroup", edge.length = max(nodeHeights(hym_sub)))
hym_sub

# Check subtree.
is.rooted(hym_sub)
is.ultrametric(hym_sub)
plot.phylo(hym_sub, cex = 0.8)
dev.off()

##--------------------------------------------------------------------------------##

#######################
### SET SIMULATIONS ###
#######################

# Set some parameters.
n_sim = 10
n_stm = 100

# Build a base Q matrix (binary).
Q <- initQ(c(0,1), c(1), diag.as = -1)

# Re-scale Q matrices (both symmetrical).
Q_sim1 <- Q*0.0001
Q_sim2 <- Q_sim1*100 # 100X

# Set character names.
char_names <- paste0("C", 1:n_sim)

# Build a list to store characters from Q1.
char_list <- setNames((vector(mode = "list", length = n_sim)), char_names)

# Sample character histories.
for (i in 1:n_sim) {
  
  # Sample histories for Q1.
  stm1 <- sim.history(hym_tree, Q_sim1, nsim = n_stm, anc = "0", message = F)
  
  for (j in 1:(n_stm)) {
    
    # Get root state of subtree.
    root_st <- names(stm1[[j]]$maps[[8]][1])
    
    # Sample histories for Q2.
    stm2 <- sim.history(hym_sub, Q_sim2, nsim = 1, anc = root_st, message = F)
    
    # Transplant clade.
    stm1[[j]]$maps[8:74] <- stm2$maps[-length(stm2$maps)]
    
  }
  
  # Store stochastic maps.
  char_list[[i]] <- stm1
  
}

# Plot some examples.
plotSimmap(char_list$C1[[1]], lwd = 3, fsize = 0.5)
plotSimmap(char_list$C2[[1]], lwd = 3, fsize = 0.5)
plotSimmap(char_list$C3[[1]], lwd = 3, fsize = 0.5)
dev.off()
