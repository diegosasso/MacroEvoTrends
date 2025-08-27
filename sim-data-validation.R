rm(list = ls())
source('R-hmm/rayDisc-multi.R')

tree <- readRDS("tree_test.RDS")
Ntip(tree)
sum(tree$edge.length)
plot.phylo(tree, cex = 0.5)
axisPhylo()
nodelabels(frame = "none", col = "blue", cex = 0.5)
#edgelabels(frame = "none", col = "blue", cex = 0.5)


#tree.paint <- paintBranches(tree, edge=167, state=1, anc.state="0")
tree.paint <- paintSubTree(tree, node=157, state=1, anc.state="0")
plot(tree.paint)

# Set some parameters.
nsim = 100
sum(tree$edge.length)*0.1

#------ Dataset 1, ONE uniform regime

Q_base <- initQ(c(0, 1), c(0.1,0.1))

# Simulate directional evolution (20 binary traits).
direct_data <- c()
for (k in 1:nsim) {
  
  # Loop until simulate an informative trait.
  x <- c()
  while (length(unique(x)) != 2) {
    anc_state <- sample(c("0", "1"), size = 1, replace = TRUE)
    #anc_state <- "0"
    x <- sim.history(tree, Q_base, anc=anc_state, nsim=1)$states
    # plot(x)
  }
  direct_data <- cbind(direct_data, x)
}

data <- cbind(rownames(direct_data), direct_data)
data_hmm <- data
data_hmm[, -1] <- ifelse(data_hmm[, -1] == "0", "0&1", "2&3")

dir <- 'R-hmm/data-validation'
file <- 'A1-ch100-reg1.rds'
path <- file.path(dir, file)
saveRDS(data_hmm, path)



#------ Dataset 2, Two regimes, one ASYM
rm(list = ls())
source('R-hmm/rayDisc-multi.R')
tree <- readRDS("tree_test.RDS")
tree.paint <- paintSubTree(tree, node=157, state=1, anc.state="0")
plot(tree.paint)
nsim = 100

Q_base <- initQ(c(0, 1), c(0.1,0.1))
Q_jump <- initQ(c(0, 1), c(0.5, 0.2))
Q_list <- list('0'=Q_base, '1'=Q_jump)
Q_list


direct_data <- c()
for (k in 1:nsim) {
  # Loop until simulate an informative trait.
  x <- c()
  while (length(unique(x)) != 2) {
    anc_state <- sample(c("0", "1"), size = 1, replace = TRUE)
    #anc_state <- "0"
    x <- sim.multiMk(tree.paint, Q = Q_list, anc = anc_state, nsim = 1)
    x_num <- as.numeric(as.character(x))
    names(x_num) <- names(x)
    
  }
  direct_data <- cbind(direct_data, x_num)
}

data <- cbind(rownames(direct_data), direct_data)
data_hmm <- data
data_hmm[, -1] <- ifelse(data_hmm[, -1] == "0", "0&1", "2&3")

dir <- 'R-hmm/data-validation'
file <- 'A2-ch100-reg2-asym.rds'
path <- file.path(dir, file)
saveRDS(data_hmm, path)


#------ Dataset 3, Two uniform regimes, ALL SYM
rm(list = ls())
source('R-hmm/rayDisc-multi.R')
tree <- readRDS("tree_test.RDS")
tree.paint <- paintSubTree(tree, node=157, state=1, anc.state="0")
plot(tree.paint)
nsim = 100

Q_base <- initQ(c(0, 1), c(0.1,0.1))
Q_jump <- initQ(c(0, 1), c(0.5, 0.5))
Q_list <- list('0'=Q_base, '1'=Q_jump)
Q_list


direct_data <- c()
for (k in 1:nsim) {
  # Loop until simulate an informative trait.
  x <- c()
  while (length(unique(x)) != 2) {
    anc_state <- sample(c("0", "1"), size = 1, replace = TRUE)
    #anc_state <- "0"
    x <- sim.multiMk(tree.paint, Q = Q_list, anc = anc_state, nsim = 1)
    x_num <- as.numeric(as.character(x))
    names(x_num) <- names(x)
    
  }
  direct_data <- cbind(direct_data, x_num)
}

data <- cbind(rownames(direct_data), direct_data)
data_hmm <- data
data_hmm[, -1] <- ifelse(data_hmm[, -1] == "0", "0&1", "2&3")

dir <- 'R-hmm/data-validation'
file <- 'A3-ch100-reg2-sym.rds'
path <- file.path(dir, file)
saveRDS(data_hmm, path)

