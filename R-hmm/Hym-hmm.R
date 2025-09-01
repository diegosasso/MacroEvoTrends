rm(list = ls())
source('R-hmm/rayDisc-multi.R')
source('R-hmm/utils-validation.R')
source('R-hmm/utils-simmap.R')
load('R-hmm/data/hym_data.RDA')


hym_mat_bin
hym_mat_recoded <- hym_mat_bin %>%
  mutate(across(-taxa, ~ recode(.x,
                                "0" = "0&1",
                                "1" = "2&3",
                                "?" = "0&1&2&3")))


Ntip(hym_tree)
plot.phylo(hym_tree, cex = 0.5)
axisPhylo()
nodelabels(frame = "none", col = "blue", cex = 0.5)
#edgelabels(frame = "none", col = "blue", cex = 0.5)



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
colnames(Q.r2.asym) <- rownames(Q.r2.asym) <- c(0:3)
Q.r2.asym_r <- Q.r2.asym[v,v]
Q.r2.asym_r
Q.r2.asym

#--- Q.r2.asym0
Q.r2.asym0 <- Q.r2.asym
Q.r2.asym0[4,2] <- NA
Q.r2.asym0[v,v]
Q.r2.asym0
#--- Q.r2.sym
Q.r2.sym <- Q.r2.asym
Q.r2.sym[4,2] <- 2
Q.r2.sym[v,v]
Q.r2.sym

#---- Q.r1.asym

Q1 <- initQ(c(0, 1), c(1, 2), diag.as = 0)
Q2 <- initQ(c('A', 'B'), c(1, 1), diag.as = 0)
Q.r1.asym <- amaSMM(Q1,Q2, diag.as = NA, non.rate.as = NA)
colnames(Q.r1.asym) <- rownames(Q.r1.asym) <- c(0:3)
Q.r1.asym
Q.r1.asym[v,v]


#---------------Validate ---------------

val.multi <- rayDISC_multi(hym_tree, hym_mat_matrix[,c(1,2)], Nchar = 1,
              rate.mat = Q.r2.asym, hmm.map = c("0&1", "2&3"),
              root.p = "maddfitz", node.states = "none",
              lewis.asc.bias = F)


val.cor <- rayDISC(hym_tree, hym_mat_matrix[,c(1,2)], rate.mat = Q.r2.asym, model='ER', root.p = "maddfitz", node.states = "none",
                           lewis.asc.bias = F)

val.multi
val.cor
#---------------Fit ---------------
hym_mat_recoded
Nchar <- 278
hym_mat_matrix <- as.matrix(hym_mat_recoded)

fit_Q.r1 <- rayDISC_multi(hym_tree, hym_mat_matrix, Nchar = Nchar,
                          rate.mat = Q.r1, hmm.map = c("0&1", "2&3"),
                          root.p = "maddfitz", node.states = "none",
                          lewis.asc.bias = TRUE)

fit_Q.r1.asym <- rayDISC_multi(hym_tree, hym_mat_matrix, Nchar = Nchar,
                          rate.mat = Q.r1.asym, hmm.map = c("0&1", "2&3"),
                          root.p = "maddfitz", node.states = "none",
                          lewis.asc.bias = TRUE)

fit_Q.r2.asym <- rayDISC_multi(hym_tree, hym_mat_matrix, Nchar = Nchar,
                               rate.mat = Q.r2.asym, hmm.map = c("0&1", "2&3"),
                               root.p = "maddfitz", node.states = "none",
                               lewis.asc.bias = TRUE)

fit_Q.r2.asym0 <- rayDISC_multi(hym_tree, hym_mat_matrix, Nchar = Nchar,
                               rate.mat = Q.r2.asym0, hmm.map = c("0&1", "2&3"),
                               root.p = "maddfitz", node.states = "none",
                               lewis.asc.bias = TRUE)

fit_Q.r2.sym <- rayDISC_multi(hym_tree, hym_mat_matrix, Nchar = Nchar,
                              rate.mat = Q.r2.sym, hmm.map = c("0&1", "2&3"),
                              root.p = "maddfitz", node.states = "none",
                              lewis.asc.bias = TRUE)

fit_Q.r1$AIC
fit_Q.r1.asym$AIC

fit_Q.r2.asym$AIC
fit_Q.r2.asym0$AIC
fit_Q.r2.sym$AIC

fit_Q.r2.sym$solution[v,v]

# madfitz
# > fit_Q.r1$AIC
# [1] 10401.83
# > fit_Q.r2.asym$AIC
# [1] 9972.841
# > fit_Q.r2.asym0$AIC
# [1] 10017.42
# > fit_Q.r2.sym$AIC
# [1] 9970.948


#  yang
# > fit_Q.r1$AIC
# [1] 10463.5
# > fit_Q.r2.asym$AIC
# [1] 10043.19
#   > fit_Q.r2.sym$AIC
# [1] 10055.81


#--flat root
# > fit_Q.r2.asym0$AIC
# [1] 10095.51
# > fit_Q.r2.asym$AIC
# [1] 10057.73
# > fit_Q.r2.asym0$AIC
# [1] 10095.51
# > fit_Q.r2.sym$AIC
# [1] 10055.81

fit_Q.r2.sym <- rayDISC_multi(hym_tree, hym_mat_matrix, Nchar = Nchar,
                              rate.mat = Q.r2.sym, hmm.map = c("0&1", "2&3"),
                              root.p = "flat", node.states = "none",
                              lewis.asc.bias = TRUE)


#--------------- Reconstruct state tips ---------------
#char.mt <- prepare_vectorized_charsSimmap(dt, tree, 200)

# arguments
char.states <- c(0:3)
Q.r2.sym
pars <- fit_Q.r2.sym$solution
Q.pars <- c(pars[1,3],
            pars[2,4],
            pars[1,2])


Q.pars
char.mt <- compute_tip_likelihoods(hym_mat_matrix[,c(1,2:200)], hym_tree, char.states=c(0:3), Q.pars=Q.pars, rate.mat=Q.r2.sym,
                                   hmm.map = c("0&1", "2&3"), root.p = "flat")

char.mt[['char1']]
ls(char.mt) %>% length()

# saveRDS(char.mt, file = 'R-hmm/data/char_mt.rds')
#--------------- Make Simmap ---------------

fit_Q.r2.sym$solution[v,v]
#tree
Nsimmap <- 100
p.root <- rep(1/4,4)

simmap_results <- simulate_multi_char_simmap(
  fit_Q = fit_Q.r2.sym,
  tree = hym_tree,
  char.mt = char.mt,
  dt = hym_mat_matrix[,c(1,2:200)],
  p.root = rep(1/4, 4),
  Nsimmap = 100
)

length(simmap_results)
simmap_results[[1]]
# saveRDS(simmap_results, file = 'R-hmm/data/simmap_results.rds')

#--------------- Calculate rates ---------------


#-----------------
# dmap<-densityMap(simmap_results[[1]], plot=T, res=100)
# 
# smap <- simmap_results[[1]][[1]]
# recode <- c('1>0', '2>b')
# smap$maps

# plot(smap)
# recode <- c('0>A', '1>A', '2>B', '3>B')
# smap.r <- recode_simmap_states(smap, recode)
# plot(smap.r)

# recode <- c('0>S', '2>S', '1>A', '3>A')
# map.recode <- recode_simmap_states(simmap_results[[1]], recode)
# dmap<-densityMap(map.recode, plot=T, res=100)
# simmap_results[1:3]

all_simmaps <- do.call("c", simmap_results)
recode <- c('0>A', '2>A', '1>S', '3>S')
all_simmaps.recode <- recode_simmap_states(all_simmaps, recode)
dmap<-densityMap(all_simmaps.recode, plot=T, res=100)
# saveRDS(dmap, file = 'R-hmm/data/dmap.rds')




#---
cols.num <- lapply(dmap$tree$maps, names) %>% unlist %>% as.numeric
hist(cols.num)
dmap$tree$maps

hist(cols.num)
cols.new <- dmap$cols
length(cols.new)
cols.new[1:1001] <- "#000000"
gradient_colors <- colorRampPalette(c("red", "blue"))(100)
#cols.new[200:250] <-"red"
#length(gradient_colors)
cols.new[120:219] <-gradient_colors
dmap.new <- setMap(dmap, cols.new)
plot(dmap.new) 



