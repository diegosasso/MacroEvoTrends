rm(list = ls())
source('R-hmm/rayDisc-multi.R')
source('R-hmm/utils-validation.R')
library(parallel)

tree <- readRDS("tree_test.RDS")
Ntip(tree)
sum(tree$edge.length)
plot.phylo(tree, cex = 0.5)
axisPhylo()
nodelabels(frame = "none", col = "blue", cex = 0.5)
#edgelabels(frame = "none", col = "blue", cex = 0.5)
#nsim=100
#Ndata=50



#--------------- MAKE Qs ---------------

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

#--------------- RUN PARALLEL ---------------


# ----- A2
dir <- 'R-hmm/data-validation'
file <- 'A2-ch100-reg2-asym.rds'
path <- file.path(dir, file)
AIC_mat <- compute_aic(path)

dir <- 'R-hmm/results'
file <- paste0('AIC_', file)
path <- file.path(dir, file)
saveRDS(AIC_mat, path)


# 
# # ----- A1
# dir <- 'R-hmm/data-validation'
# file <- 'A1-ch100-reg1.rds'
# path <- file.path(dir, file)
# AIC_mat <- compute_aic(path)
# 
# dir <- 'R-hmm/results'
# file <- 'AIC_A1-ch100-reg1.rds'
# path <- file.path(dir, file)
# saveRDS(AIC_mat, path)
# 
# 
# 
# 
# # ----- A3
# dir <- 'R-hmm/data-validation'
# file <- 'A3-ch100-reg2-sym.rds'
# path <- file.path(dir, file)
# AIC_mat <- compute_aic(path)
# 
# dir <- 'R-hmm/results'
# file <- paste0('AIC_', file)
# path <- file.path(dir, file)
# saveRDS(AIC_mat, path)

# #--------------- RUN PARALLEL ---------------
# #--------------- SETTINGS -------------------
# dir <- "R-hmm/data-validation"
# files <- list.files(dir, pattern = "A.*\\.rds$", full.names = TRUE)
# files
# ncores <- detectCores() - 1  # use all but one core
# 
# #--------------- SETTINGS -------------------
# 
# ncores <- detectCores() - 1 
# cl <- makeCluster(ncores)
# clusterExport(cl, varlist = c("tree", "Q.r1", "Q.r2.asym", "Q.r2.sym",
#                               "rayDISC_multi"))  # export needed vars
# # Make all helper functions available to each worker
# clusterEvalQ(cl, {
#   source("R-hmm/utils-hmm.R")
#   source("R-hmm/rayDISC-multi.R")
# })
# 
# 
# #A1.AIC <- compute_aic(path, Q.r1, Q.r2.asym, Q.r2.sym)
# AIC_list <- parLapply(cl, files, compute_aic)
# stopCluster(cl)
# 
# # Combine into one big matrix
# AIC_all <- do.call(rbind, AIC_list)
# 
# # Save results
# #saveRDS(AIC_all, file = output_file)
# 
# 
# 
# #------------------ Fit A1 with reg1
# dir <- 'R-hmm/data-validation'
# file <- 'A1-ch100-reg1.rds'
# path <- file.path(dir, file)
# data_hmm <- readRDS(path)
# 
# A1.AIC <- compute_aic(path, Q.r1, Q.r2.asym, Q.r2.sym)
# 
# #----
# fit_Q.r1 <- rayDISC_multi(tree, data_hmm[[i]], Nchar = nsim, rate.mat= Q.r1, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
# fit_Q.r2.asym <- rayDISC_multi(tree, data_hmm[[i]], Nchar = nsim, rate.mat= Q.r2.asym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
# fit_Q.r2.sym <- rayDISC_multi(tree, data_hmm[[i]], Nchar = nsim, rate.mat= Q.r2.sym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
# 
# 
# AIC <- c(fit_Q.r1$AIC,
#          fit_Q.r2.asym$AIC,
#          fit_Q.r2.sym$AIC)
# 
# 
# 
# 
# dAIC.r2 <- fit_Q.r1$AIC-fit_Q.r2.asym$AIC
# dAIC.r2.asym <-  fit_Q.r1$AIC-fit_Q.r2.sym$AIC
# 
# 
# fit_Q.r1 # best
# fit_Q.r2.asym 
# fit_Q.r2.sym
# 
# #------------------ Fit A2 with reg2 ASYM
# dir <- 'R-hmm/data-validation'
# file <- 'A2-ch100-reg2-asym.rds'
# path <- file.path(dir, file)
# data_hmm <- readRDS(path)
# 
# fit_Q.r1 <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r1, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
# fit_Q.r2.asym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.asym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
# fit_Q.r2.sym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.sym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
# 
# fit_Q.r1 
# fit_Q.r2.asym # best
# fit_Q.r2.sym
# 
# 
# #------------------ Fit A3 with reg2 SYM
# dir <- 'R-hmm/data-validation'
# file <- 'A3-ch100-reg2-sym.rds'
# path <- file.path(dir, file)
# data_hmm <- readRDS(path)
# 
# fit_Q.r1 <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r1, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
# fit_Q.r2.asym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.asym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
# fit_Q.r2.sym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Q.r2.sym, hmm.map = c('0&1', '2&3'), root.p='flat', node.states="none", lewis.asc.bias = TRUE)
# 
# fit_Q.r1 
# fit_Q.r2.asym # best
# fit_Q.r2.sym # best
# 

