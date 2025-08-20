
# Load packages.
library(ape)
library(phytools)
library(corHMM)
library(rphenoscate)
library(ontophylo)

# Import some functions.
source('R/utils.R')

# Import tree.
tree <- readRDS("tree_test.RDS")
#tree <- readRDS("hym_tree.RDS")

# Simulate tree.
#tree <- pbtree(n = 100, scale = 1, b = 1, d = 0)
plot.phylo(tree, cex = 0.5)
nodelabels(frame = "none", col = "red", cex = 0.5)
edgelabels(frame = "none", col = "blue", cex = 0.5)
dev.off()

# Save tree.
#saveRDS(tree, "tree_test.RDS")

# Simulate directional selection.
# Set focal edge.
focal_edge <- 5

# Set base Q-matrix.
Q_base <- initQ(c(0, 1), c(0.01,0.01))
Q_base

# Set Q-matrix of shift.
# Note that Q_jump might be reordered depending on initial state.
Q_jump <- initQ(c(0, 1), c(10, 1e-10))
Q_jump

# Set some parameters.
nsim = 20

# Simulate directional evolution (20 binary traits).
direct_data <- c()
for (k in 1:nsim) {
  
  # Loop until simulate an informative trait.
  x <- c()
  while (length(unique(x)) != 2) {
    
    x <- simDirectional(tree, focal_edge, Q_base, Q_jump, nsim = 1)
    
  }
  
  direct_data <- cbind(direct_data, x)
  
}

# Set trait labels.
colnames(direct_data) <- paste0("C",1:nsim)
direct_data

# Check tip data for traits.
plot.phylo(tree, show.tip.label = F)
nodelabels(frame = "none", col = "red", cex = 0.5)
edgelabels(frame = "none", col = "blue", cex = 0.5)
tiplabels(tip = (direct_data[,1] == "1"), pch = 19, col = "purple", cex = 0.5)
dev.off()

# Set some parameters.
nstm = 100

# Sample stochastic histories for each simulated trait.
# List of fitted models.
fitted_models <- setNames(vector(mode = "list", length = nsim), paste0("C",1:nsim))

# List of stochastic maps.
smaps <- setNames(vector(mode = "list", length = nsim), paste0("C",1:nsim))

i = 1
# Loop over all traits.
for (i in 1:nsim){
  
  cat(paste0("\n", "Working on: ", "C", i, ": ", Sys.time(), "\n"))
  
  # Set character vector.
  char <- cbind(rownames(direct_data), direct_data[,i])
  
  # Build set of models.
  mods <- c("ER", "SYM", "ARD")
  
  # Fit models.
  fit_Q <- setNames(vector(mode = "list", length = length(mods)), mods)
  
  j = 1
  for (j in 1:length(mods)) {
    
    fit_Q[[j]] <- corHMM(phy = tree, data = char, model = mods[[j]], rate.cat = 1, root.p = "yang")
    
  }
  
  # Get best model.
  w <- geiger::aicw(sapply(fit_Q, function(x) x$AICc))[,3]
  
  # Get Q matrix.
  Q <- fit_Q[[min(which(w == max(w)))]]$solution
  
  # Store fitted Q matrix.
  fitted_models[[i]] <- Q
  
  # Sample stochastic maps.
  smaps[[i]] <- makeSimmap(tree = tree, data = char, model = Q, rate.cat = 1, nSim = nstm)
  
}

# Check some reconstructions.
densityMap(smaps$C1)
densityMap(smaps$C2)
densityMap(smaps$C3)
densityMap(smaps$C4)
dev.off()

# Set some parameters.
res = 1000

# Discretize all tress.
stm_discr <- lapply(smaps, function(x) discr_Simmap_all(x, res = res) )

# Amalgamate characters.
stm_amalg <- paramo.list(names(smaps), tree.list = stm_discr, ntrees = nstm)

# Merge state categories across branches.
stm_merg <- merge_tree_cat_list(stm_amalg)

# Discretize reference tree.
tree_discr <- discr_Simmap(tree, res = res)

# Calculate Hamming distances.
cat(paste0("\n", "Starting calculating hamming distances: ", Sys.time(), "\n"))
path_hm <- path_hamming_over_trees_KDE(stm_merg)
cat(paste0("\n", "Finished calculating hamming distances: ", Sys.time(), "\n"))

# Make path data.
path_data <- make_data_NHPP_KDE_Markov_kernel(path_hm, add.psd = F)

# Check for branches with no changes.
sapply(path_data, length)

# Estimate bandwidth.
bdw <- estimate_band_W(tree_discr, path_data, band.width = "bw.nrd")
bdw <- mean(bdw)

# Kernel Density Estimator (KDE).
edge_KDE <- estimate_edge_KDE(tree_discr, Path.data = path_data, h = bdw)

# Calculate smoothing and normalize KDE data.
edge_KDE$Maps.mean.loess <- loess_smoothing_KDE(tree_discr, edge_KDE)
edge_KDE$Maps.mean.loess.norm <- normalize_KDE(tree_discr, edge_KDE$Maps.mean.loess)

# Calculate the lambda statistics.
lambda_post <- posterior_lambda_KDE(stm_merg)
lambda_post

# Get the posterior distribution.
edge_KDE$lambda.mean <- make_postPois_KDE(edge_KDE$Maps.mean.norm, lambda_post, lambda.post.stat = "Mean")
edge_KDE$lambda.mean.loess <- make_postPois_KDE(edge_KDE$Maps.mean.loess.norm, lambda_post, lambda.post.stat = "Mean")

# Make data for contmaps.
nhpp_lambda_mean <- make_contMap_KDE(tree_discr, edge_KDE$lambda.mean.loess)

# Make data for edge profiles.
edge_profs_lambda_mean <- edge_profiles4plotting(tree_discr, edge_KDE$lambda.mean.loess)

# Plot. need to fix plotting function!
# edgeplot(nhpp_lambda_mean, edge_profs_lambda_mean)

# Load package.
library(tidyverse)
library(ggplot2)

# Make plot.
png("edgeplot_direct.png", units = "in", width = 7, height = 7, res = 300)

# Get tree height.
Tmax <- max(nodeHeights(nhpp_lambda_mean$tree))

# Set plot layout.
layout(matrix(c(1,2),ncol = 1), heights = c(2,1))

# Plot contmap.
plot.contMap(nhpp_lambda_mean, lwd = 3, outline = F, legend = F, ftype = "off",
             plot = F, mar = c(0.1, 3.45, 0.1, 0.35))

# Plot edgeplot.
plot_edgeprof <-
  
  ggplot(data = edge_profs_lambda_mean, aes(x = X-Tmax, y =  Y, group = edge.id, color = Y)) +
  
  geom_line(alpha = 1, linewidth = 0.5) + 
  
  scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = 0.7)) ) +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        plot.margin = unit(c(2.3,0.87,0.1,0.1), 'cm'),
        legend.position = 'none') +
  
  xlab('time') + ylab('rate') +
  
  scale_x_continuous(limits = c(-round(Tmax, 5), 0),
                     breaks = -1*seq(from = 0, to = Tmax, by = Tmax/5) %>% round(5),
                     labels = seq(from = 0, to = Tmax, by = Tmax/5) %>% round(5) ) +
  scale_y_continuous(limits = c(-0.5, round(max(edge_profs_lambda_mean$Y)*1.2, 3))) +
  coord_cartesian(expand = FALSE)

vp <- grid::viewport(height = unit(0.5,"npc"), width = unit(1, "npc"), just = c("left",'top'), y = 0.5, x = 0)

print(plot_edgeprof, vp = vp)

title(main = paste0("Branch rates"), font.main = 2, line = -0.5, cex.main = 0.8)

dev.off()

# Save data.
#saveRDS(stm, "smaps.RDS")
#saveRDS(stm_amalg, "stm_amalg.RDS")
#save(direct_data, fitted_models, path_hm, path_data, edge_KDE, file = "test.RDA")

# Get a morphospace.
MD <- suppressWarnings(MultiScale.simmap(stm_merg[[1]]))

# Plot. need to fix plotting function!
# mds_plot(MD)

# Load package.
library(tidyverse)
library(ggplot2)

# Add tip ids.
MD$Points <- mutate(MD$Points, tip.id = c(1:nrow(MD$Points)))

# Get tree height.
Tmax <- max(MD$Points$time)

# Get time slice.
Tslice = max(MD$Points$time)

# Filter data according to time slice.
MD$Points <- filter(MD$Points, MD$Points$time <= Tslice)

# Generate MDS plot.
GG <-
  
  ggplot(MD$Points) +
  
  geom_segment(data = MD$Lines, aes(x = MD$Lines$start.V1, y = MD$Lines$start.V2, xend = MD$Lines$end.V1, yend = MD$Lines$end.V2),
               colour = "red", linewidth = 0.3, linetype = 1, alpha = 0.3) +
  
  geom_point(aes(x = MD$Points$V1, y = MD$Points$V2, color = MD$Points$time, group = MD$Points$tip.id), alpha = 0.8, size = 1.6) +
  
  scale_colour_gradient(low = "purple",  high = "green", limits = c(0, Tmax),
                        breaks = seq(0, Tmax, Tmax/5), labels = round(seq(0, Tmax, Tmax/5),2) %>% rev) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),) +
  
  xlab("Coor1") + ylab("Coor2") +
  guides(color = guide_colourbar(title = "Time Myr", reverse = F, draw.ulim = T, draw.llim = T))


# Make plot.
png("morphospace_direct.png", units = "in", width = 7, height = 7, res = 300)
GG
dev.off()
