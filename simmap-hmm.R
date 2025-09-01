rm(list = ls())
source('R-hmm/rayDisc-multi.R')
source('R-hmm/utils-validation.R')
source('R-hmm/utils-simmap.R')

tree <- readRDS("tree_test.RDS")
Ntip(tree)
sum(tree$edge.length)
plot.phylo(tree, cex = 0.5)
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

#--- Q.r2.sym
Q.r2.sym <- Q.r2.asym
Q.r2.sym[4,2] <- 2
Q.r2.sym[v,v]
Q.r2.sym

#---- chech results
dir <- 'R-hmm/results'
file <- 'AIC_A2-ch100-reg2-asym.rds'
path <- file.path(dir, file)
AIC_A2 <- readRDS(path)
AIC_A2 <- as_tibble(AIC_A2)
AIC_A2

hist(AIC_A2$Q.r2.asym - AIC_A2$Q.r2.sym)
sum((AIC_A2$Q.r2.asym - AIC_A2$Q.r2.sym) < -2)
which((AIC_A2$Q.r2.asym - AIC_A2$Q.r2.sym) < -10)
# we take the 1st

#--------------- Fit A2 ---------------
dir <- 'R-hmm/data-validation'
file <- 'A2-ch100-reg2-asym.rds'
path <- file.path(dir, file)
data_hmm <- readRDS(path)
dt <- data_hmm[[1]] 

fit_Q.r2.asym <- rayDISC_multi(tree, dt, Nchar = ncol(dt)-1,
                               rate.mat = Q.r2.asym, hmm.map = c("0&1", "2&3"),
                               root.p = "flat", node.states = "none",
                               lewis.asc.bias = TRUE)

fit_Q.r2.asym
AIC_A2

#--------------- Reconstruct state tips ---------------
#char.mt <- prepare_vectorized_charsSimmap(dt, tree, 200)

# arguments
char.states <- c(0:3)
pars <- fit_Q.r2.asym$solution
Q.pars <- c(pars[1,3],
pars[2,4],
pars[4,3],
pars[4,2])


char.mt <- compute_tip_likelihoods(dt, tree, char.states=c(0:3), Q.pars=Q.pars, rate.mat=Q.r2.asym,
                                   hmm.map = c("0&1", "2&3"), root.p = "flat")

char.mt[['char1']]
ls(char.mt) %>% length()

#--------------- Make Simmap ---------------

fit_Q.r2.asym$solution[v,v]
tree
Nsimmap <- 100
p.root <- rep(1/4,4)

simmap_results <- simulate_multi_char_simmap(
  fit_Q = fit_Q.r2.asym,
  tree = tree,
  char.mt = char.mt,
  dt = dt,
  p.root = rep(1/4, 4),
  Nsimmap = 100
)

length(simmap_results)
simmap_results[[1]]

#--------------- Calculate rates ---------------

direc.changes <- c('0>1', '2>3')
#get_per_branch_rate(simmap_results[[1]], direc.changes)
r.dir <- sum_per_branch_rates(simmap_results, direc.changes)
hist(r.dir)

indirec.changes <- c('1>0', '3>2')
r.indir <- sum_per_branch_rates(simmap_results, indirec.changes)
hist(r.indir)

plot(r.dir, r.indir)

# Create a color scale (blue → red)
colors <- colorRampPalette(c("blue", "yellow", "red"))(100)
# Map r.dir values to color indices
r.norm <- (r.dir- min(r.dir)) / (max(r.dir) - min(r.dir))  # Normalize [0,1]
branch_colors <- colors[ceiling(r.norm * 99) + 1]
# Plot tree colored by r.dir
plot(tree, show.tip.label = FALSE, edge.color = branch_colors, edge.width = 2)

#-----
direc.changes <- c('1>3', '3>1')
#get_per_branch_rate(simmap_results[[1]], direc.changes)
r.dir <- sum_per_branch_rates(simmap_results, direc.changes)
hist(r.dir)

indirec.changes <- c('0>2', '2>0')
r.indir <- sum_per_branch_rates(simmap_results, indirec.changes)
hist(r.indir)

plot(r.dir, r.indir)

hist(r.dir-r.indir)
sum((r.dir-r.indir) > 0)

# Create a color scale (blue → red)
colors <- colorRampPalette(c("blue", "yellow", "red"))(100)
# Map r.dir values to color indices
r.norm <- (r.dir- min(r.dir)) / (max(r.dir) - min(r.dir))  # Normalize [0,1]
branch_colors <- colors[ceiling(r.norm * 99) + 1]
# Plot tree colored by r.dir
plot(tree, show.tip.label = FALSE, edge.color = branch_colors, edge.width = 2)

#---

#------ plot

smap <- simmap_results[[1]][[1]]
plot(smap)
smap$maps
smap$edge.length
r.dir
max(r.dir)
which(r.dir>15)

#direc.changes <- c('0>1', '2>3', '1>3', '3>1')
direc.changes <- c( '0>1', '2>3')
#get_per_branch_rate(simmap_results[[1]], direc.changes)
r.dir <- sum_per_branch_rates(simmap_results, direc.changes)
hist(r.dir)

#indirec.changes <- c('1>0', '3>2', '0>2', '2>0')
indirec.changes <- c('1>0', '3>2')
#get_per_branch_rate(simmap_results[[1]], direc.changes)
r.indir <- sum_per_branch_rates(simmap_results, indirec.changes)
hist(r.indir)

hist(r.dir-r.indir)
r.diff <- (r.dir-r.indir)
r.diff[r.diff<0] <- -10

# Create a color scale (blue → red)
colors <- colorRampPalette(c("blue", "yellow", "red"))(100)
# Map r.dir values to color indices
r.norm <- (r.diff- min(r.diff)) / (max(r.diff) - min(r.diff))  # Normalize [0,1]
branch_colors <- colors[ceiling(r.norm * 99) + 1]
# Plot tree colored by r.dir
plot(tree, show.tip.label = FALSE, edge.color = branch_colors, edge.width = 2)



#--- indirect
# Create a color scale (blue → red)
colors <- colorRampPalette(c("blue", "yellow", "red"))(100)
# Map r.dir values to color indices
r.norm <- (r.indir - min(r.indir)) / (max(r.indir) - min(r.indir))  # Normalize [0,1]
branch_colors <- colors[ceiling(r.norm * 99) + 1]
# Plot tree colored by r.dir
plot(tree, show.tip.label = FALSE, edge.color = branch_colors, edge.width = 2)

plot.phylo(tree, cex = 0.5)
edgelabels(frame = "none", col = "blue", cex = 0.5)


#-----------------
dmap<-densityMap(simmap_results[[1]], plot=T, res=100)

smap <- simmap_results[[1]][[1]]
recode <- c('1>0', '2>b')
smap$maps

plot(smap)
recode <- c('0>A', '1>A', '2>B', '3>B')
smap.r <- recode_simmap_states(smap, recode)
plot(smap.r)

recode <- c('0>S', '2>S', '1>A', '3>A')
map.recode <- recode_simmap_states(simmap_results[[1]], recode)
dmap<-densityMap(map.recode, plot=T, res=100)
simmap_results[1:3]

all_simmaps <- do.call("c", simmap_results)
recode <- c('0>S', '2>S', '1>A', '3>A')
all_simmaps.recode <- recode_simmap_states(all_simmaps, recode)
dmap<-densityMap(all_simmaps.recode, plot=T, res=100)
dmap$states
dmap$tree$maps
cols.num <- lapply(dmap$tree$maps, names) %>% unlist %>% as.numeric
hist(cols.num)
dmap$cols
dmap1 <- dmap
cols.new <- dmap1$cols
# Create the red → blue gradient for positions 1 to 500
gradient_colors <- colorRampPalette(c("red", "blue"))(31)
plot(gradient_colors)
# Create black for positions 501 to 1001
black_colors <- rep("#000000", 1001-600)

# Combine into a single palette
custom_palette <- c(gradient_colors, black_colors)

# Optional: check the first and last few colors
head(custom_palette)   # red end
tail(custom_palette)   # all black
# length(cols.new)
cols.new[1:1001] <- "#000000"
#length(cols.new[540:570])
cols.new[540:570] <- gradient_colors
dmap1$cols <- cols.new
#dmap1$cols[1:1001] <-custom_palette
plot(dmap1)
dmap.new <- setMap(dmap, cols.new)
plot(dmap.new)
