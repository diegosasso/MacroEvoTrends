library(expm)
library(ape)
library(phytools)
library(corHMM)
library(rphenoscate)
library(ontophylo)
library(dplyr)
source('R/lik-hamming.R')


tree <- readRDS("tree_test.RDS")
# plot.phylo(tree, cex = 0.5)
# nodelabels(frame = "none", col = "blue", cex = 0.5)
# edgelabels(frame = "none", col = "blue", cex = 0.5)


tree.paint <- paintBranches(tree, edge=167, state=1, anc.state="0")
plot(tree.paint)

Q_base <- initQ(c(0, 1), c(0.01,0.01))
Q_base
Q_jump <- initQ(c(0, 1), c(10, 1e-10))
Q_jump
Q_list <- list('0'=Q_base, '1'=Q_jump)

char <- sim.multiMk(tree.paint, Q_list, anc='0', nsim=2)
char[] <- lapply(char, function(x) as.numeric(as.character(x)))
head(char)
str(char)
#saveRDS(char, 'RDS/char.rds')

Q.est <- initQ(c(0, 1), c(1, 1), diag.as = NA)
taxa <- cbind(rownames(char), char[,1])
recon <- rayDISC(tree, taxa, rate.mat= Q.est, model='ARD', node.states="none", lewis.asc.bias = TRUE)
recon

#--------------------------------------
source('R-hmm/rayDisc-multi.R')

tree <- readRDS("tree_test.RDS")
char <- readRDS('RDS/char.rds')
data <- cbind(rownames(char), char)
str(data)
Q.est <- initQ(c(0, 1), c(1, 1), diag.as = NA)

fit_multi <- rayDISC_multi(tree, data, Nchar = 2, rate.mat= Q.est, p=c(0.1), root.p='flat', node.states="none", lewis.asc.bias = F)

taxa1 <- cbind(rownames(char), char[,1])
recon1 <- rayDISC(tree, taxa1,  p=c(0.1), rate.mat= Q.est, root.p=c(0.5,0.5), model='ER', node.states="none", lewis.asc.bias = F)
taxa2 <- cbind(rownames(char), char[,2])
recon2 <- rayDISC(tree, taxa2,  p=c(0.1), rate.mat= Q.est, root.p=c(0.5,0.5), model='ER', node.states="none", lewis.asc.bias = F)
recon1$loglik + recon2$loglik
fit_multi$loglik

# --- Assym
Q.est2 <- initQ(c(0, 1), c(1, 2), diag.as = NA)

fit_multi <- rayDISC_multi(tree, data, Nchar = 2, rate.mat = Q.est2,  p=c(0.1, 0.01), root.p='flat', node.states="none", lewis.asc.bias = F)

taxa1 <- cbind(rownames(char), char[,1])
recon1 <- rayDISC(tree, taxa1, rate.mat = Q.est2, p=c(0.1,0.01), root.p=c(0.5,0.5), model='ARD', node.states="none", lewis.asc.bias = F)
taxa2 <- cbind(rownames(char), char[,2])
recon2 <- rayDISC(tree, taxa2, rate.mat = Q.est2, p=c(0.1,0.01), root.p=c(0.5,0.5), model='ARD', node.states="none", lewis.asc.bias = F)
recon1$loglik + recon2$loglik
fit_multi$loglik

#---- Optimize
Q.est <- initQ(c(0, 1), c(1, 1), diag.as = NA)
fit_multi <- rayDISC_multi(tree, data, Nchar = 2, rate.mat= Q.est, root.p='flat', node.states="none", lewis.asc.bias = F)


Q.est2 <- initQ(c(0, 1), c(1, 2), diag.as = NA)
fit_multi <- rayDISC_multi(tree, data, Nchar = 2, rate.mat = Q.est2,  root.p='flat', node.states="none", lewis.asc.bias = F)

recon1 <- rayDISC(tree, taxa1, rate.mat = Q.est2,  root.p=c(0.5,0.5), model='ARD', node.states="none", lewis.asc.bias = F)


#---- Optimize many chars
tree <- readRDS("tree_test.RDS")
Ntip(tree)
# plot.phylo(tree, cex = 0.5)
# nodelabels(frame = "none", col = "blue", cex = 0.5)
# edgelabels(frame = "none", col = "blue", cex = 0.5)


tree.paint <- paintBranches(tree, edge=167, state=1, anc.state="0")
plot(tree.paint)

Q_base <- initQ(c(0, 1), c(0.02,0.02))
Q_base
Q_jump <- initQ(c(0, 1), c(0.5, 1e-10))
Q_jump
Q_list <- list('0'=Q_base, '1'=Q_jump)
Q_list <- list('0'=Q_base, '1'=Q_base)
Q_list
# Set some parameters.
nsim = 50

# Simulate directional evolution (20 binary traits).
direct_data <- c()
for (k in 1:nsim) {
  
  # Loop until simulate an informative trait.
  x <- c()
  while (length(unique(x)) != 2) {
    #anc_state <- sample(c("0", "1"), size = 1, replace = TRUE)
    anc_state <- "0"
    x <- sim.multiMk(tree.paint, Q = Q_list, anc = anc_state, nsim = 1)
    x_num <- as.numeric(as.character(x))
    names(x_num) <- names(x)
    
  }
  
  direct_data <- cbind(direct_data, x_num)
  
}

data <- cbind(rownames(direct_data), direct_data)

Q.est <- initQ(c(0, 1), c(1, 1), diag.as = NA)
fit_multi1 <- rayDISC_multi(tree, data, Nchar = nsim, rate.mat= Q.est, root.p=c(1,0), node.states="none", lewis.asc.bias = F)
fit_multi1_L <- rayDISC_multi(tree, data, Nchar = nsim, rate.mat= Q.est, root.p=c(1,0), node.states="none", lewis.asc.bias = TRUE)

Q.est2 <- initQ(c(0, 1), c(1, 2), diag.as = NA)
fit_multi2 <- rayDISC_multi(tree, data, Nchar = nsim, rate.mat= Q.est2, root.p=c(1,0), node.states="none", lewis.asc.bias = F)
fit_multi2_L <- rayDISC_multi(tree, data, Nchar = nsim, rate.mat= Q.est2, root.p=c(1,0), node.states="none", lewis.asc.bias = TRUE)

fit_multi1_L
fit_multi2_L

fit_multi1
fit_multi2




 
