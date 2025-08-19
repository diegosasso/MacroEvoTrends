library(expm)
library(ape)
library(phytools)
library(rphenoscate)
source('R/utils.R')


# Simulate tree
#tree<-pbtree(n=5, scale=1, b=1, d=0)
tree.txt <- "(t1:1,(t2:0.752918711,(t3:0.03401850044,(t4:0.01232941458,t5:0.01232941458):0.02168908586):0.7189002105):0.247081289);"
tree <- read.tree(text=tree.txt)
plot(tree)
nodelabels()
edgelabels()

#--- Simulate directional selection, TOY TREE
focal_edge <- 4
Q_base <- initQ(c(0, 1), c(1,1))
# not that Q_jump might be reodered depending on init state
Q_jump <- initQ(c(0, 1), c(10, 1e-10))
Q_jump

#simDirectional_one_char(tree, focal_edge, Q_base, Q_jump)
simDirectional(tree, focal_edge, Q_base, Q_jump, nsim=20)


#--- Simulate directional selection, Bigger tree

tree<-pbtree(n=100, scale=1, b=1, d=0)
plot(tree)
edgelabels()

# e.g., 2
focal_edge <- 1
Q_base <- initQ(c(0, 1), c(.01,.01))
Q_jump <- initQ(c(0, 1), c(10, 1e-10))
Q_base
Q_jump

chars <- simDirectional(tree, focal_edge, Q_base, Q_jump, nsim=20)
chars[,2]

