library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)
library(viridis)
source("R/utils-ST.R")
source('R/region_class.R')

stm_merg <-readRDS("data/stm_merg.RDS")
tree <-readRDS("data/hym_tree.RDS")
plot(tree)

# branch rate
rates_list <- lapply(stm_merg, get_branch_rate_across_smaps)
rates <- do.call(cbind, rates_list)
head(rates)
rate_total <- apply(rates, 2, sum)
rate_total


# Min-max scale each column to [0, 1]
rates_scaled <- apply(rates, 2, function(col) {
  (col - min(col)) / (max(col) - min(col))
})
head(rates_scaled)
max(rates_scaled[,1])


#----- Overrepresentation
enr <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")

enr.branch <- apply(enr, 1, sum)
enr.boolean <- (enr.branch > 0)*1
hist(enr.branch)

plot_branch_colored_tree(tree, enr.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')


#------ Convert to binary
enr.bin <- apply(enr, 1, function(x) paste0(x, collapse = ''))
enr.bin
unique(enr.bin)
sum(enr.bin > 0)
enr.boolean

# region_class is sourced from the file
region_class
enr_class <- enrich_next_level(enr, region_class)
enr_class <- enr_class[c("head", "mesosoma", "metasoma", "legs",  "wings")]
head(enr_class)
head(enr)

# Get enrichment for 5 regions
bin.enr <- apply(enr_class, 1, function(x) paste0(x, collapse = ''))
bin.enr[grepl("^0+$", bin.enr)] <- ""
unique(bin.enr) %>% length()
table(bin.enr)
sum(enr_class > 0)


#------------- Plot numerical

plot_branch_colored_tree(tree, enr.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')


edgelabels(bin.enr, frame = "none", col = "red", cex = 0.6)
#edgelabels( frame = "none", col = "red", cex = 0.6)


#------------- Plot navajo

plot_branch_colored_tree(tree, enr.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')


edge_coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)

for (i in seq_along(bin.enr)) {
  state <- bin.enr[i]
  if (state == "") next
  
  # edge endpoints
  x0 <- edge_coords$xx[tree$edge[i, 1]]
  y0 <- edge_coords$yy[tree$edge[i, 1]]
  x1 <- edge_coords$xx[tree$edge[i, 2]]
  y1 <- edge_coords$yy[tree$edge[i, 2]]
  
  # edge midpoint (x only), y = child’s y (keeps glyphs aligned)
  xm <- (x0 + x1) / 2
  ym <- y1
  
  # binary string → vector
  bits <- as.integer(strsplit(state, "")[[1]])
  n <- length(bits)
  
  # rectangle dimensions
  box_w <- max(edge_coords$xx) * 0.01
  box_h <- max(edge_coords$yy) * 0.01
  gap   <- box_w * 0.1
  
  total_w <- n * (box_w + gap)
  x_start <- xm - total_w/2
  
  # draw boxes
  for (j in seq_len(n)) {
    xleft   <- x_start + (j-1) * (box_w + gap)
    xright  <- xleft + box_w
    ybottom <- ym - box_h/2
    ytop    <- ym + box_h/2
    
    rect(xleft, ybottom, xright, ytop,
         col = ifelse(bits[j] == 1, "red", "white"),
         border = "black", lwd = 0.5)
  }
}



#------------ Underrepresented

under <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "under", return_type = "binary")

under.branch <- apply(under, 1, sum)
under.branch
under.boolean <- (under.branch > 0)*1
hist(under.branch)

plot_branch_colored_tree(tree, under.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')

#------ Convert to binary
under.bin <- apply(under, 1, function(x) paste0(x, collapse = ''))
under.bin
unique(under.bin)
sum(under.boolean > 0)
under.boolean

# region_class is sourced from the file
region_class
under_class <- enrich_next_level(under, region_class)
under_class <- under_class[c("head", "mesosoma", "metasoma", "legs",  "wings")]
head(under_class)
head(under)

# Get under for 5 regions
under.enr <- apply(under_class, 1, function(x) paste0(x, collapse = ''))
under.enr[grepl("^0+$",under.enr)] <- ""
unique(under.enr) %>% length()
table(under.enr)
sum(under.enr > 0)

#------------- Plot numerical

plot_branch_colored_tree(tree, under.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')


edgelabels(under.enr, frame = "none", col = "red", cex = 0.6)
#edgelabels( frame = "none", col = "red", cex = 0.6)


#------------- Plot navajo

plot_branch_colored_tree(tree, under.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')


edge_coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)

for (i in seq_along(under.enr)) {
  state <- under.enr[i]
  if (state == "") next
  
  # edge endpoints
  x0 <- edge_coords$xx[tree$edge[i, 1]]
  y0 <- edge_coords$yy[tree$edge[i, 1]]
  x1 <- edge_coords$xx[tree$edge[i, 2]]
  y1 <- edge_coords$yy[tree$edge[i, 2]]
  
  # edge midpoint (x only), y = child’s y (keeps glyphs aligned)
  xm <- (x0 + x1) / 2
  ym <- y1
  
  # binary string → vector
  bits <- as.integer(strsplit(state, "")[[1]])
  n <- length(bits)
  
  # rectangle dimensions
  box_w <- max(edge_coords$xx) * 0.01
  box_h <- max(edge_coords$yy) * 0.01
  gap   <- box_w * 0.1
  
  total_w <- n * (box_w + gap)
  x_start <- xm - total_w/2
  
  # draw boxes
  for (j in seq_len(n)) {
    xleft   <- x_start + (j-1) * (box_w + gap)
    xright  <- xleft + box_w
    ybottom <- ym - box_h/2
    ytop    <- ym + box_h/2
    
    rect(xleft, ybottom, xright, ytop,
         col = ifelse(bits[j] == 1, "red", "white"),
         border = "black", lwd = 0.5)
  }
}

#------------ Tip representation

plot_branch_colored_tree(tree, enr.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')


#edgelabels(bin.enr, frame = "none", col = "red", cex = 0.6)
edgelabels( frame = "none", col = "red", cex = 0.6)


enr.masked <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "masked")

br <- enr.masked[,1]

tip_index <- 1
#tree$tip.label[77]
tip_path <- get_tip_path(tree, tip_index)

br[tip_path]

#------
enr.branch <- apply(enr, 1, sum)
enr.boolean <- (enr.branch > 0)*1
hist(enr.branch)

plot_branch_colored_tree(tree, enr.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')