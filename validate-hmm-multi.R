rm(list = ls())
source('R-hmm/rayDisc-multi.R')

tree <- readRDS("tree_test.RDS")
Ntip(tree)
sum(tree$edge.length)
plot.phylo(tree, cex = 0.5)
axisPhylo()
nodelabels(frame = "none", col = "blue", cex = 0.5)
#edgelabels(frame = "none", col = "blue", cex = 0.5)
nsim=100

#------------------ Make Qs
Q1 <- initQ(c(0, 1), c(1, 1), diag.as = 0)
Q2 <- initQ(c('A', 'B'), c(1, 1), diag.as = 0)
Qsmm <- amaSMM(Q1,Q2, diag.as = NA, non.rate.as = NA)

#--- Q.r1
Q.r1 <- Qsmm
#colnames(Q.r1) <- rownames(Q.r1) <- c(0:3)

#--- Q.r2.asym
Q.r2.asym <- Qsmm
Q.r2.asym[2,4] <- 2
Q.r2.asym[4,2] <- 4
Q.r2.asym[1,2] <- 3
Q.r2.asym[2,1] <- 3
Q.r2.asym[3,4] <- 3
Q.r2.asym[4,3] <- 3
v <- c(1,3,2,4)
Q.r2.asym_r <- Q.r2.asym[v,v]
Q.r2.asym_r
Q.r2.asym

#--- Q.r2.sym
Q.r2.sym <- Q.r2.asym
Q.r2.sym[4,2] <- 2
Q.r2.sym[v,v]
Q.r2.sym



#------------------ Fit A1 with reg1
dir <- 'R-hmm/data-validation'
file <- 'A1-ch100-reg1.rds'
path <- file.path(dir, file)
data_hmm <- readRDS(path)


fit_Q.r1 <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r1, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
fit_Q.r2.asym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.asym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
fit_Q.r2.sym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.sym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)

fit_Q.r1 # best
fit_Q.r2.asym 
fit_Q.r2.sym

#------------------ Fit A2 with reg2 ASYM
dir <- 'R-hmm/data-validation'
file <- 'A2-ch100-reg2-asym.rds'
path <- file.path(dir, file)
data_hmm <- readRDS(path)

fit_Q.r1 <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r1, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
fit_Q.r2.asym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.asym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
fit_Q.r2.sym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.sym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)

fit_Q.r1 
fit_Q.r2.asym 
fit_Q.r2.sym


#------------------ Fit A3 with reg2 SYM
dir <- 'R-hmm/data-validation'
file <- 'A3-ch100-reg2-sym.rds'
path <- file.path(dir, file)
data_hmm <- readRDS(path)

fit_Q.r1 <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r1, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
fit_Q.r2.asym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.asym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
fit_Q.r2.sym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.sym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)

fit_Q.r1 
fit_Q.r2.asym 
fit_Q.r2.sym