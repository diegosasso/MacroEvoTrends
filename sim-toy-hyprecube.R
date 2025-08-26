
# Load packages.
library(ape)
library(phytools)
library(corHMM)
library(rphenoscate)
library(ontophylo)
library(dplyr)

# Import some functions.
source('R/utils.R')

#load("test.RDA")

# Import tree.
tree <- readRDS("tree_test.RDS")
#tree <- readRDS("hym_tree.RDS")
plot(tree)

# Simtrees# Simulate tree.
#tree <- pbtree(n = 100, scale = 1, b = 1, d = 0)
plot.phylo(tree, cex = 0.5)
#nodelabels(frame = "none", col = "red", cex = 0.5)
edgelabels(frame = "none", col = "blue", cex = 0.5)
dev.off()

# Save tree.
#saveRDS(tree, "tree_test.RDS")

# Simulate directional selection.
# Set focal edge.
focal_edge <- 154

# Set base Q-matrix.
Q_base <- initQ(c(0, 1), c(0.01,0.01))
Q_base

# Set Q-matrix of shift.
# Note that Q_jump might be reordered depending on initial state.
Q_jump <- initQ(c(0, 1), c(0.1, 1e-10))
Q_jump

# Set some parameters.
nsim = 50

# Simulate directional evolution (20 binary traits).
direct_data <- c()
for (k in 1:nsim) {
  
  # Loop until simulate an informative trait.
  x <- c()
  while (length(unique(x)) != 2) {
    
    x <- simDirectional(tree, focal_edge, Q_base, Q_jump, nsim = 1)
    
  }
  
  direct_data <- cbind(direct_data, x)
  
}

# Set trait labels.
colnames(direct_data) <- paste0("C",1:nsim)
direct_data



# Set some parameters.
nstm = 100

# Sample stochastic histories for each simulated trait.
# List of fitted models.
fitted_models <- setNames(vector(mode = "list", length = nsim), paste0("C",1:nsim))

# List of stochastic maps.
smaps <- setNames(vector(mode = "list", length = nsim), paste0("C",1:nsim))

i = 1
# Loop over all traits.
for (i in 1:nsim){
  
  cat(paste0("\n", "Working on: ", "C", i, ": ", Sys.time(), "\n"))
  
  # Set character vector.
  char <- cbind(rownames(direct_data), direct_data[,i])
  
  # Build set of models.
  mods <- c("SYM", "ARD")
  
  # Fit models.
  fit_Q <- setNames(vector(mode = "list", length = length(mods)), mods)
  
  j = 1
  for (j in 1:length(mods)) {
    
    fit_Q[[j]] <- corHMM(phy = tree, data = char, model = mods[[j]], rate.cat = 1, root.p = "yang")
    
  }
  
  # Get best model.
  w <- geiger::aicw(sapply(fit_Q, function(x) x$AICc))[,3]
  
  # Get Q matrix.
  Q <- fit_Q[[min(which(w == max(w)))]]$solution
  
  # Store fitted Q matrix.
  fitted_models[[i]] <- Q
  
  # Sample stochastic maps.
  smaps[[i]] <- makeSimmap(tree = tree, data = char, model = Q, rate.cat = 1, nSim = nstm)
  
}



# Set some parameters.
res = 1000

# Discretize all tress.
stm_discr <- lapply(smaps, function(x) discr_Simmap_all(x, res = res) )

# Amalgamate characters.
stm_amalg <- paramo.list(names(smaps), tree.list = stm_discr, ntrees = nstm)

# Merge state categories across branches.
stm_merg <- merge_tree_cat_list(stm_amalg)

# Discretize reference tree.
tree_discr <- discr_Simmap(tree, res = res)

lapply(stm_merg, function(x) x$maps[[5]] )

saveRDS(stm_merg, 'RDS/stm_merg_toy-ST.RDS')
