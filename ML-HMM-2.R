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

#------ Asym Data
Q_base <- initQ(c(0, 1), c(0.05,0.05))
Q_jump <- initQ(c(0, 1), c(0.4, 0.1))
Q_list <- list('0'=Q_base, '1'=Q_jump)
Q_list

# Set some parameters.
nsim = 200

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
data_hmm <- data
data_hmm[, -1] <- ifelse(data_hmm[, -1] == "0", "0&1", "2&3")


#------ Sym Data
Q_base <- initQ(c(0, 1), c(0.05,0.05))
Q_jump <- initQ(c(0, 1), c(0.4, 0.4))
Q_list <- list('0'=Q_base, '1'=Q_jump)
Q_list

# Set some parameters.
nsim = 200

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

data_sym <- cbind(rownames(direct_data), direct_data)
data_hmm_sym <- data_sym
data_hmm_sym[, -1] <- ifelse(data_hmm_sym[, -1] == "0", "0&1", "2&3")

#------- Assym Fit
Q1 <- initQ(c(0, 1), c(1, 1), diag.as = 0)
Q2 <- initQ(c('A', 'B'), c(1, 1), diag.as = 0)
Qsmm <- amaSMM(Q1,Q2, diag.as = NA, non.rate.as = NA)
Qsmm[2,4] <- 2
Qsmm[4,2] <- 4
Qsmm[1,2] <- 3
Qsmm[2,1] <- NA
Qsmm[3,4] <- 3
Qsmm[4,3] <- NA
v <- c(1,3,2,4)
Qsmm_r <- Qsmm[v,v]
Qsmm_r
Qsmm
colnames(Qsmm) <- rownames(Qsmm) <- c(0:3)
Qsmm

source('R-hmm/rayDisc-multi.R')
fit_asym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Qsmm, hmm.map = c('0&1', '2&3'), root.p=c(1,0,0,0), node.states="none", lewis.asc.bias = TRUE)
fit_sym <- rayDISC_multi(tree, data_hmm_sym, Nchar = nsim, rate.mat= Qsmm, hmm.map = c('0&1', '2&3'), root.p=c(1,0,0,0), node.states="none", lewis.asc.bias = TRUE)

fit_asym
fit_sym
fit_asym$solution[v,v]
fit_sym$solution[v,v]
Q_list

# simmap <- makeSimmap(tree=tree, data=data[,c(1,2)], model=fit_asym$solution, root.p="yang", rate.cat=2, nSim=1)
# 
# data(primates)
# phy <- primates[[1]]
# phy <- multi2di(phy)
# data <- primates[[2]]
# 
# ##run corhmm
# MK <- corHMM(phy, data, 1)
# ancRECON
# ##get simmap from corhmm solution
# model <- MK$solution
# simmap <- makeSimmap(tree=phy, data=data, model=model, rate.cat=1, nSim=1, nCores=1)
# plotSimmap(simmap[[1]])
# 
# corHMM:::GetTipStateBruteForce
#------------------------- 2 rates
Qsmm2 <- Qsmm
Qsmm2[4,2] <- 2
Qsmm2[v,v]
fit_r2_asym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Qsmm2, hmm.map = c('0&1', '2&3'), root.p=c(1,0,0,0), node.states="none", lewis.asc.bias = TRUE)
fit_r2_sym <- rayDISC_multi(tree, data_hmm_sym, Nchar = nsim, rate.mat= Qsmm2, hmm.map = c('0&1', '2&3'), root.p=c(1,0,0,0), node.states="none", lewis.asc.bias = TRUE)

fit_r2_asym
fit_r2_sym
fit_r2_asym$solution[v,v]
fit_r2_sym$solution[v,v]


#------------------------- Sym Fit
Qsmm1 <- amaSMM(Q1,Q2, diag.as = NA, non.rate.as = NA)
colnames(Qsmm1) <- rownames(Qsmm1) <- c(0:3)
Qsmm1
fit_r1_asym <- rayDISC_multi(tree, data_hmm, Nchar = nsim, rate.mat= Qsmm1, hmm.map = c('0&1', '2&3'), root.p=c(1,0,0,0), node.states="none", lewis.asc.bias = TRUE)
fit_r1_sym <- rayDISC_multi(tree, data_hmm_sym, Nchar = nsim, rate.mat= Qsmm1, hmm.map = c('0&1', '2&3'), root.p=c(1,0,0,0), node.states="none", lewis.asc.bias = TRUE)

fit_r1_asym
fit_r1_sym
fit_r1_asym$solution[v,v]
fit_r1_sym$solution[v,v]

#------
fit_asym
fit_r2_asym
fit_r1_asym

fit_sym
fit_r2_sym
fit_r1_sym
