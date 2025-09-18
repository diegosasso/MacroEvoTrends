# Load libraries.
library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)
source("R/utils-ST.R")


stm_merg <-readRDS("data/stm_merg.RDS")
tree <-readRDS("data/hym_tree.RDS")
plot(tree)

# branch rate
rates_list <- lapply(stm_merg, get_branch_rate_across_smaps)
rates <- do.call(cbind, rates_list)
head(rates)
rate_total <- apply(rates, 2, sum)
rate_total

rate_total_branch <- apply(rates, 1, sum)


# Min-max scale each column to [0, 1]
rates_scaled <- apply(rates, 2, function(col) {
  (col - min(col)) / (max(col) - min(col))
})
head(rates_scaled)
max(rates_scaled[,1])

# ---- 5 second level BD  ----
stm_merg2 <-  readRDS("data/stm_merg2.RDS")
stm_merg2
rates_list2 <- lapply(stm_merg2, get_branch_rate_across_smaps)
rates2 <- do.call(cbind, rates_list2)
colMeans(rates2 == 0)
names(stm_merg2)

rates2_scaled <- apply(rates2, 2, function(col) {
  (col - min(col)) / (max(col) - min(col))
})
head(rates2_scaled)

library(viridis)

enr <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")
enr2 <- get_enrichment_mask(rates2_scaled, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")
dis <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "under", return_type = "binary")

# Alll
all.enr <- apply(enr, 1, sum)
all.dis <- apply(dis, 1, sum)
all.enr
all.dis

hist(all.enr)

plot_branch_colored_tree(tree, (all.enr > 0)*1 , palette = viridis(15), 
                         title = "all", legend_title = 'Cols')

plot_branch_colored_tree(tree, (all.enr<4 & all.enr>0)*1 , palette = viridis(15), 
                         title = "all", legend_title = 'Cols')

plot_branch_colored_tree(tree, (all.enr>4)*1 , palette = viridis(15), 
                         title = "all", legend_title = 'Cols')

two.regs <- (all.enr>0) + (all.dis>0)
two.regs <-(two.regs==2)*1

plot_branch_colored_tree(tree, two.regs , palette = viridis(15), 
                         title = "two.regs", legend_title = 'Cols')

plot_branch_colored_tree(tree, (all.dis>0)*1, palette = viridis(15), 
                         title = "Dis", legend_title = 'Cols')


r1 <- (all.enr > 0)*1
#r2 <- (two.regs>0)*1
r3 <- (all.dis>0)*1
rr <- paste0(r1,  r3, sep='')
unique(rr)
mapping <- c("00" = 1, "01" = 2, "10" = 3, "11" = 4)
rr_num <- unname(mapping[rr])

plot_branch_colored_tree(tree, rr_num, palette =c('black', 'green', 'red', 'blue'), 
                         title = "4 regimes, 15bd", legend_title = 'Cols')

# Head
head <- apply(enr[,c(1:2)], 1, sum)

cbind(enr[,c(1:2)], enr2[,1])

plot_branch_colored_tree(tree, head , palette = viridis(15), 
                         title = "Head", legend_title = 'Cols')

# Thorax
thorax <- apply(enr[,c(3:8)], 1, sum)
table(thorax)
plot_branch_colored_tree(tree, thorax , palette = viridis(15), 
                         title = "Head", legend_title = 'Cols')

#----- stm2
enr2 <- get_enrichment_mask(rates2_scaled, lower_q = 0.1, upper_q = 0.9, direction = "over", return_type = "binary")
dis2 <- get_enrichment_mask(rates2_scaled, lower_q = 0.05, upper_q = 0.95, direction = "under", return_type = "binary")

all.enr2 <- apply(enr2, 1, sum)
#------
bin.enr <- apply(enr, 1, function(x) paste0(x, collapse = ''))
bin.enr
unique(bin.enr)
sum(all.enr > 0)
all.enr

region_class
enr

# --- Example input ---
# enr is your enrichment matrix (branches × body regions)
# region_class from before (with region → group)

# Convert to data.frame for easier handling
enr_df <- as.data.frame(enr)
enr_df$branch <- seq_len(nrow(enr_df))  # keep branch ID

# Reshape to long format
library(reshape2)
enr_long <- melt(enr_df, id.vars = "branch",
                 variable.name = "region",
                 value.name = "enriched")

source('adj-matrix.R')
# Join with region_class to map each region to a group
enr_long <- merge(enr_long, region_class, by = "region")

# Collapse: for each branch × group, if any enriched==1 → mark 1
enr_collapsed <- aggregate(enriched ~ branch + group, data = enr_long,
                           FUN = function(x) as.integer(any(x > 0)))

# Reshape back to wide matrix: branches × groups
enr_class <- reshape(enr_collapsed,
                     idvar = "branch",
                     timevar = "group",
                     direction = "wide")

# Clean column names
colnames(enr_class) <- sub("enriched\\.", "", colnames(enr_class))

# Remove branch col if not needed
rownames(enr_class) <- enr_class$branch
enr_class <- enr_class[ , -1]

enr_class
head(enr_class)
head(enr)

bin.enr <- apply(enr_class, 1, function(x) paste0(x, collapse = ''))
unique(bin.enr)
table(bin.enr)
sum(enr_class > 0)

plot(tree)
axisPhylo()
edgelabels(bin.enr, frame = "none", col = "red", cex = 0.6)
#-------------

bin.enr[grepl("^0+$", bin.enr)] <- ""

plot_branch_colored_tree(tree, (all.enr > 0)*1 , palette = viridis(15), 
                         title = "all", legend_title = 'Cols')

edgelabels( frame = "none", col = "red", cex = 0.6)
edgelabels(bin.enr, frame = "none", col = "red", cex = 0.6)

#_-----

plot_branch_colored_tree(tree, (all.enr > 0)*1 , palette = viridis(15), 
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

# # Example tree
# set.seed(1)
# tree <- rtree(5)

# Get edge heights (start and end times for each edge)
H <- nodeHeights(tree)  # matrix with (start, end) for each edge
max(nodeHeights(tree))
max(H[,2])
# Say you want to test for time point t
t <-0.1

# Find edges overlapping t
#edge_ids <- which(H[,1] >= t & H[,2] <= t)
edge_ids <- which(t >= H[,1] & t<=H[,2] )
# Edge IDs (rows of tree$edge correspond to H)
edge_ids
tree$edge[edge_ids, ]

#--------------------------

library(phytools)
library(reshape2)
library(ggplot2)

enr <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")

# Example usage:
denom = c("active", "enriched")
plot_enrichment_stacked(tree, enr, n_points = 1000,  denom = "active")
plot_enrichment_stacked(tree, enr, n_points = 1000,  denom = "enriched")
colnames(enr)
p1 <- plot_enrichment_stacked(tree, enr[,c(1:2)], n_points = 10000)
p2 <-plot_enrichment_stacked(tree, enr[,c(3:8)], n_points = 10000)
p3 <-plot_enrichment_stacked(tree, enr[,c(9:10)], n_points = 1000)
p4 <-plot_enrichment_stacked(tree, enr[,c(11:13)], n_points = 1000)
p5 <-plot_enrichment_stacked(tree, enr[,c(14:15)], n_points = 1000)

p1 <- plot_enrichment_stacked(tree, enr[,c(1:2)], n_points = 10000, denom = "enriched")
p2 <-plot_enrichment_stacked(tree, enr[,c(3:8)], n_points = 10000, denom = "enriched")
p3 <-plot_enrichment_stacked(tree, enr[,c(9:10)], n_points = 1000, denom = "enriched")
p4 <-plot_enrichment_stacked(tree, enr[,c(11:13)], n_points = 1000, denom = "enriched")
p5 <-plot_enrichment_stacked(tree, enr[,c(14:15)], n_points = 1000, denom = "enriched")


dis <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "under", return_type = "binary")

p1 <- plot_enrichment_stacked(tree, dis[,c(1:2)], n_points = 1000, denom = "enriched")
p2 <-plot_enrichment_stacked(tree, dis[,c(3:8)], n_points = 1000, denom = "enriched")
p3 <-plot_enrichment_stacked(tree, dis[,c(9:10)], n_points = 100, denom = "enriched")
p4 <-plot_enrichment_stacked(tree, dis[,c(11:13)], n_points = 100, denom = "enriched")
p5 <-plot_enrichment_stacked(tree, dis[,c(14:15)], n_points = 100, denom = "enriched")

library(cowplot)

plot_grid(
  p1, p2, p3, p4, p5,
  ncol = 1,  # one column → stacked vertically
  align = "v"
)


#-------- Entropy measure diversity enrcihed BRs (i.e. the adaptive diversification)
enr <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")

p1 <- plot_enrichment_entropy(tree, enr, type = "Neff", n_points = 10000)
tail(p1$data)
plot_enrichment_entropy(tree, enr, type = "Neff_scaled", n_points = 10000)
#plot_enrichment_entropy(tree, enr, type = "weighted", n_points = 10000)

plot_enrichment_entropy(tree, enr[,8, drop=F], type = "Neff_scaled", n_points = 10000)

plot_enrichment_entropy(tree, enr[,c(1:2)], type = "Neff_scaled", n_points = 10000)
plot_enrichment_entropy(tree, enr[,c(3:8)], type = "Neff_scaled", n_points = 10000)
plot_enrichment_entropy(tree, enr[,c(9:10)], type = "Neff_scaled", n_points = 10000)
plot_enrichment_entropy(tree, enr[,c(11:13)], type = "Neff_scaled", n_points = 10000)
plot_enrichment_entropy(tree, enr[,c(14:15)], type = "Neff_scaled", n_points = 10000)

plot(p1$data$entropy, p2$data$entropy, typ='l')
plot(p3$data$entropy, p2$data$entropy, typ='l')

#-- Pyrate
pyr <- read.csv('data/pyrate.csv', header = F)

plot(pyr$V1, pyr$V2, type='l')

x <- plot_enrichment_entropy(tree, enr, type = "Neff_scaled", n_points = 27)$data


plot(x$time, x$entropy, type = 'l')
plot(pyr$V1, pyr$V2, type='l')

en <- x$entropy
di <- pyr$V2
plot(di, en)
di[1] <- di[2]
plot(pyr$V1, di, type='l')
log(di)
di <- di - min(di) + 1e-6   # add offset to make min slightly > 0
log(di)
log(en)
en[1] <- 1e-6
log(en)

plot(di, en)
plot(log(di), log(en))
plot(log(di), en)

plot(di, en)
abline(lm(en ~ di), col = "red", lwd = 2)
summary(lm(en ~ di))
cor.test(en, di, method = "pearson")

plot(log(di), en)
abline(lm(en ~ log(di)), col = "red", lwd = 2)
cor.test(en, log(di), method = "pearson")


x <- plot_enrichment_entropy(tree, enr, type = "Neff_scaled", n_points = 10000)$data
di_time <- pyr$V1
di
di_time

library(ggplot2)

# dataframe for entropy curve
df_entropy <- data.frame(
  time = x$time,
  entropy = x$entropy
)

# dataframe for di
df_di <- data.frame(
  time = di_time,
  di = di
)

ggplot() +
  # entropy curve (dense line)
  geom_line(data = df_entropy, aes(x = time, y = entropy, color = "Neff_scaled"), size = 1) +
  
  # di curve (line + points at sparse time steps)
  geom_line(data = df_di, aes(x = time, y = di / max(di, na.rm = TRUE) * max(df_entropy$entropy, na.rm = TRUE),
                              color = "di"), size = 1) +
  geom_point(data = df_di, aes(x = time, y = di / max(di, na.rm = TRUE) * max(df_entropy$entropy, na.rm = TRUE),
                               color = "di"), size = 2) +
  
  scale_y_continuous(
    name = "Neff_scaled",
    sec.axis = sec_axis(~ . * max(di, na.rm = TRUE) / max(df_entropy$entropy, na.rm = TRUE),
                        name = "di")
  ) +
  scale_color_manual(values = c("Neff_scaled" = "blue", "di" = "red")) +
  scale_x_reverse() +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.title.y.left  = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "red")
  ) +
  labs(x = "Time from root (Ma)", color = "Metric")


#-------- Enrichment strength mesuares simultenous strong mode of evolution across BRs

enr <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")
dis <- get_enrichment_mask(rates_scaled, lower_q = 0.05, upper_q = 0.95, direction = "under", return_type = "binary")

plot_enrichment_strength(tree, enr, type = "scaled", n_points = 1000)
plot_enrichment_strength(tree, enr, type = "not_scaled", n_points = 1000)
tail(p1$data)

plot_enrichment_entropy(tree, enr, type = "Neff_scaled", n_points = 1000)
plot_enrichment_entropy(tree, dis, type = "Neff_scaled", n_points = 1000)
plot_enrichment_strength(tree, dis, type = "scaled", n_points = 1000)

library(ggplot2)


# ---- Collect the three metrics ----

# Neff_scaled for enr
df1 <- plot_enrichment_entropy(tree, enr, type = "Neff_scaled", n_points = 1000 )$data
names(df1)[names(df1) == "entropy"] <- "metric_value"
df1$metric <- "Neff_scaled_enr"

# Strength for enr
df2 <- plot_enrichment_strength(tree, enr, type = "not_scaled", n_points = 1000 )$data
names(df2)[names(df2) == "strength"] <- "metric_value"
df2$metric <- "Strength_scaled_enr"

# Neff_scaled for dis
df3 <- plot_enrichment_entropy(tree, dis, type = "Neff_scaled", n_points = 1000)$data
names(df3)[names(df3) == "entropy"] <- "metric_value"
df3$metric <- "Neff_scaled_dis"

# Merge
df_all <- rbind(df1, df2, df3)

# ---- Scaling for strength metric ----
range1 <- range(c(df1$metric_value, df3$metric_value), na.rm = TRUE)  # Neff ranges
range2 <- range(df2$metric_value, na.rm = TRUE)                      # Strength ranges
scale_factor <- diff(range1) / diff(range2)

df_all$metric_scaled <- with(df_all,
                             ifelse(metric == "Strength_scaled_enr",
                                    metric_value * scale_factor,
                                    metric_value)
)

# ---- Plot ----
# p <- ggplot(df_all, aes(x = time, y = metric_scaled, color = metric)) +
#   geom_line(na.rm = TRUE) +
#   theme_minimal() +
#   scale_x_reverse() +
#   scale_y_continuous(
#     name = "Neff_scaled",
#     sec.axis = sec_axis(~ . / scale_factor,
#                         name = "Enrichment_strength_scaled")
#   ) +
#   labs(x = "Time from root") +
#   theme(legend.position = "bottom")

p <- ggplot(df_all, aes(x = time, y = metric_scaled, color = metric)) +
  geom_line(size = 1) +                                # thicker lines
  theme_classic(base_size = 14) +                      # clean Nature-like theme
  scale_x_reverse(expand = c(0, 0)) +                  # root left, no padding
  scale_y_continuous(
    name = "Neff_scaled",
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Enrichment_strength_scaled")
  ) +
  labs(x = "Time from root") +
  scale_color_manual(
    values = c("Neff_scaled_dis" = "#1f77b4",          # blue
               "Neff_scaled_enr" = "#2ca02c",          # green
               "Strength_scaled_enr" = "#d62728")      # red
  ) +
  theme(
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "top",                           # move legend up
    legend.box = "horizontal",
    panel.grid = element_blank()                       # remove gridlines
  )

print(p)




p <- ggplot(df_all, aes(x = time, y = metric_scaled, color = metric, fill = metric)) +
  geom_line(size = 1) +
  geom_area(alpha = 0.2, position = "identity") +
  theme_classic(base_size = 14) +
  scale_x_reverse(
    expand = c(0, 0),                              # no gap on x
    breaks = seq(0, max(df_all$time, na.rm = TRUE), by = 50)  # ticks every 50
  ) +
  scale_y_continuous(
    expand = c(0, 0),                              # no gap on y
    name = "Neff_scaled",
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Enrichment_strength_scaled")
  ) +
  labs(x = "Time from root") +
  scale_color_manual(
    values = c("Neff_scaled_dis" = "#1f78b4",
               "Neff_scaled_enr" = "#33a02c",
               "Strength_scaled_enr" = "#e31a1c")
  ) +
  scale_fill_manual(
    values = c("Neff_scaled_dis" = "#1f78b4",
               "Neff_scaled_enr" = "#33a02c",
               "Strength_scaled_enr" = "#e31a1c")
  ) +
  theme(
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    panel.grid = element_blank()
  )

print(p)


library(deeptime)

p <- ggplot(df_all, aes(x = time, y = metric_scaled, color = metric, fill = metric)) +
  geom_line(size = 1) +
  geom_area(alpha = 0.2, position = "identity") +
  theme_classic(base_size = 14) +
  scale_x_reverse(
    expand = c(0, 0),
    breaks = seq(0, max(df_all$time, na.rm = TRUE), by = 50)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    name = "Neff_scaled",
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Enrichment_strength_scaled")
  ) +
  labs(x = "Time (Ma)") +
  scale_color_manual(
    values = c("Neff_scaled_dis" = "#1f78b4",
               "Neff_scaled_enr" = "#33a02c",
               "Strength_scaled_enr" = "#e31a1c")
  ) +
  scale_fill_manual(
    values = c("Neff_scaled_dis" = "#1f78b4",
               "Neff_scaled_enr" = "#33a02c",
               "Strength_scaled_enr" = "#e31a1c")
  ) +
  theme(
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    panel.grid = element_blank()
  ) +
  coord_geo(
    dat = "periods",     # can also use "epochs" or "stages"
    pos = "bottom",      # bottom only
    size = 3             # label size
  )

print(p)

