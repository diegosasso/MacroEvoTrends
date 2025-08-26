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
corHMM:::dev.raydisc

