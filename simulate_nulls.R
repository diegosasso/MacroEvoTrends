
# Load packages.
library(ape)
library(phytools)
library(ontophylo)

# Import data from 'empirical' analyses.
load("test.RDA")

# Import tree.
tree <- readRDS("tree_test.RDS")
plot(tree)
# Set some parameters.
nsim = 20
nstm = 100
res = 500

# List of stochastic maps.
smaps <- setNames(vector(mode = "list", length = nsim), paste0("C",1:nsim))

# Loop over all traits.
i=1
for (i in 1:nsim) {
  
  cat(paste0("\n", "Working on: ", "C", i, ": ", Sys.time(), "\n"))
  
  # Build Q matrix from models fitted to the observed data.
  Q <- fitted_models[[i]]
  diag(Q) <- -rowSums(Q, na.rm = T)
  Q[is.na(Q)] <- 0
  colnames(Q) <- rownames(Q) <- 0:(dim(Q)[1] - 1)
  
  # Simulate character histories.
  smaps[[i]] <- sim.history(tree, Q = Q, nsim = nstm)
  
}

# Discretize all tress.
stm_discr <- lapply(smaps, function(x) discr_Simmap_all(x, res = res) )

# Amalgamate characters.
stm_amalg <- paramo.list(names(smaps), tree.list = stm_discr, ntrees = nstm)

# Merge state categories across branches.
stm_merg <- merge_tree_cat_list(stm_amalg)

# Rename objects.
smaps_null <- smaps
stm_amalg_null <- stm_amalg

# Save data.
saveRDS(smaps_null, "stm_null.RDS")
saveRDS(stm_amalg_null, "stm_amalg_null.RDS")
