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

# tot states per BR
load('data/char_num.RDA')
tot_states <- lapply(anat_ent_state, function(x) sum(get_len(x)) ) %>% unlist
tot_states


rates_norm <- apply(rates, 1, function(x) x/tot_states) %>% t
head(rates_norm)
head(rates)

# Min-max scale each column to [0, 1]
# rates_scaled <- apply(rates, 2, function(col) {
#   (col - min(col)) / (max(col) - min(col))
# })
# head(rates_scaled)
# max(rates_scaled[,1])
# hist(rates_scaled[,1])
# hist(rates[,1])

#----- Overrepresentation
enr <- get_enrichment_mask(rates, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")

enr.branch <- apply(enr, 1, sum)
#enr.branch1 <- enr.branch
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

#------------ Tip overrepresentation

plot_branch_colored_tree(tree, enr.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')


#edgelabels(bin.enr, frame = "none", col = "red", cex = 0.6)
edgelabels( frame = "none", col = "red", cex = 0.6)


enr.masked <- get_enrichment_mask(rates_norm, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "masked")
#enr <- get_enrichment_mask(rates, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")
#rbind(enr.masked[,1],enr[,1])

tip.enr <- get_tip_enrichment(tree, enr.masked)
tip.enr[tip.enr>0] %>% hist()
max(tip.enr)
min(tip.enr)

# fix zeros
eps <- min(tip.enr[tip.enr > 0]) / 10  # smaller than smallest positive
tip.enr <- log(tip.enr + eps)
head(tip.enr)
max(tip.enr)
min(tip.enr)

hist(tip.enr)
tip.enr[tip.enr>log(eps)] %>% hist()
#-----------------
library(RColorBrewer)

# # custom Nature-like palette
# nature_pal <- colorRampPalette(c("#4575b4", "#fee090", "#d73027"))
# 
# # preview
# image(matrix(seq(-1, 1, length.out = 100), ncol = 1),
#       col = nature_pal(100),
#       axes = FALSE, main = "Nature-style sequential palette")
# 
# val_min <- min(tip.enr, na.rm = TRUE)
# val_max <- max(tip.enr, na.rm = TRUE)
# cols_fun <- nature_pal(100)
# 
# map_to_color <- function(vals) {
#   breaks <- seq(val_min, val_max, length.out = 101)
#   cols_fun[cut(vals, breaks = breaks, include.lowest = TRUE)]
# }
#-----
make_custom_palette <- function(val_min, val_max, midpoint,
                                low_col = "#4575b4", 
                                mid_col = "#fee090", 
                                high_col = "#d73027",
                                n = 200) {
  # ensure even n
  if (n %% 2 == 1) n <- n + 1
  
  lower_breaks <- seq(val_min, midpoint, length.out = n/2)
  upper_breaks <- seq(midpoint, val_max, length.out = n/2)
  
  # remove duplicate midpoint
  breaks <- c(lower_breaks, upper_breaks[-1])
  
  pal <- colorRampPalette(c(low_col, mid_col, high_col), bias=.5)(length(breaks)-1)
  
  list(breaks = breaks, pal = pal, midpoint = midpoint)
}


map_to_color <- function(vals, pal_info) {
  cols <- pal_info$pal[cut(vals, breaks = pal_info$breaks, include.lowest = TRUE)]
  #cols[vals == log(eps)] <- "#ffffff"  # force exact zeros to white
  #cols[vals == 0] <- "#ffffff"
  #cols[vals == 0] <- "lightgrey"
  cols[vals == 0] <- "aliceblue"  #"skyblue1" "#deebf7" 
  cols
}

#val_min <- tip.enr[tip.enr>log(eps)] %>% min #min(tip.enr, na.rm = TRUE)
val_min <- min(tip.enr[tip.enr > 0]) #min(tip.enr, na.rm = TRUE)
val_max <- max(tip.enr, na.rm = TRUE)

format(6e-5, scientific = FALSE)
midpoint <- 4e-5  # you choose this



pal_info <- make_custom_palette(val_min, val_max, midpoint, n = 64)
pal_info$pal <- gplots::rich.colors(64, palette="temperature", alpha=1.0, rgb=FALSE, plot=F)
pal_info$pal[1:8] <- pal_info$pal[9]
bampal<- BAMMtools::assignColorBreaks(c(tip.enr[tip.enr > 0]), spex=NULL, NCOLORS = 64, method="jenks")
unique(bampal) %>% length()
pal_info$breaks <- unique(bampal)
unique(pal_info$breaks) %>% length()
cols <- map_to_color(tip.enr, pal_info)
map_to_color(tip.enr[1,1], pal_info)

format(12e-5, scientific = FALSE)
format(2e-5, scientific = FALSE)
vals <- seq(val_min , 12e-5, length.out=100)

cols <- map_to_color(vals, pal_info)
# simple strip plot
par(mar = c(2, 4, 2, 2))
plot(vals, rep(1, length(vals)), 
     col = cols, pch = 15, cex = 3,
     xlab = "Value", ylab = "", yaxt = "n")

which(tip.enr>-6)
which(test>-6)
map_to_color(tip.enr[110], pal_info)
#---------
library(ape)



# ---- Plot tree ----
#branch_colors <- ifelse(enr.boolean == 1, "#FF3300FF", "black")
branch_colors <- ifelse(enr.boolean == 1, "#2ca02c", "black")
branch_widths <- ifelse(enr.boolean == 1, 2.5, 1.3)  # thicker if enriched
# "purple3" "forestgreen""deepskyblue3"

plot(tree, edge.color = branch_colors, show.tip.label = TRUE, 
     cex = 0.5, label.offset = 90, edge.width = branch_widths)

lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Tip coordinates
tip_y <- lastPP$yy[1:Ntip(tree)]
tip_x <- max(lastPP$xx)

# ---- Square sizes ----
dy <- median(diff(sort(tip_y)))
box_size_y <- dy * 1.2   # adjust multiplier for bigger/smaller squares

# Scaling factor: keep squares true
pin <- par("pin")  # plot dimensions (inches)
xrange <- diff(lastPP$x.lim)
yrange <- diff(lastPP$y.lim)
x_per_inch <- xrange / pin[1]
y_per_inch <- yrange / pin[2]
box_size_x <- box_size_y * (x_per_inch / y_per_inch)

# ---- Gap between tree and heatmap ----
gap <- 0.005 * max(lastPP$xx)

# ---- Define groups of columns ----
# Example: 3 groups 1:2, 3:8, 9:15
groups <- list(1:2, 3:8, 9:10, 11:13, 14:15)

# Compute group offsets
gap_size <- box_size_x * 0.4   # adjust for wider/narrower gaps
group_offsets <- rep(0, ncol(tip.enr))
for (g in seq_along(groups)) {
  if (g > 1) {
    group_offsets[groups[[g]]] <- (g-1) * gap_size
  }
}

# ---- Draw glyphs ----
for (i in 1:Ntip(tree)) {
  values <- as.numeric(tip.enr[i, ])
  n <- length(values)
  for (j in 1:n) {
    val <- values[j]
    col <- map_to_color(val, pal_info)
    
    # Add group offset to x position
    xleft   <- tip_x + gap + (j-1) * box_size_x + group_offsets[j]
    xright  <- xleft + box_size_x
    ybottom <- tip_y[i] - box_size_y/2
    ytop    <- tip_y[i] + box_size_y/2
    
    rect(xleft, ybottom, xright, ytop,
         col = col, border = "black", lwd = 0.4)
  }
}



#------------- Plot navajo
#plot(tree, edge.color = branch_colors, show.tip.label = TRUE, cex = 0.5, label.offset = 90, edge.width = branch_widths)

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
  # box_w <- max(edge_coords$xx) * 0.01
  # box_h <- max(edge_coords$yy) * 0.01
  box_w <- max(edge_coords$xx) * 0.011
  box_h <- max(edge_coords$yy) * 0.011
  gap   <- box_w * 0
  
  total_w <- n * (box_w + gap)
  x_start <- xm - total_w/2
  
  # draw boxes
  for (j in seq_len(n)) {
    xleft   <- x_start + (j-1) * (box_w + gap)
    xright  <- xleft + box_w
    ybottom <- ym - box_h/2
    ytop    <- ym + box_h/2
    
    rect(xleft, ybottom, xright, ytop,
         col = ifelse(bits[j] == 1, "seagreen3", "white"),
         border = "black", lwd = 0.5)
  }
}

"mediumorchid3"
•	"darkorange2"
•	"goldenrod2"
•	"gold"
