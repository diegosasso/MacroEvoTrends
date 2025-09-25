library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)
source("R/utils-ST.R")


# ====================== Read Data ====================== #

# ---- 15 body regions (first order) ----
stm_merg <- readRDS("data/stm_merg.RDS")
# branch rate
rates_list <- lapply(stm_merg, get_branch_rate_across_smaps)
rates <- do.call(cbind, rates_list)
head(rates)
# proportion of zeros
colMeans(rates == 0)

# total rate per body regions (is it correct ????)
rate_total <- apply(rates, 2, sum) !!!!! treelength !!!! / nrow(rates) # should it be devided by N branches
rate_total
# ---- Number of states per BR ----
load('data/char_num.RDA')
#anat_ent_state
# log of states per character per BR
log_states <- lapply(anat_ent_state, function(x) sum(log(get_len(x))) ) %>% unlist
log_states
tot_states <- lapply(anat_ent_state, function(x) sum(get_len(x)) ) %>% unlist
tot_states

rate_total_norm <- rate_total/tot_states
# rate_total_norm <- exp(log(rate_total)-log_states)
rate_total_norm

plot(rate_total, tot_states)
plot(rate_total, log_states)
points(rate_total, tot_states, col='red')

cor.test(rate_total, log_states, method = "pearson")
cor.test(rate_total, tot_states, method = "pearson")


# ====================== Pairwise correlation: 15 body regions (first order) ====================== #

#. !!!! Quantiles should be trimmed AFTER LOG!!!!!!!!!!!
corr.pearson.log <-run_all_pairwise_correlations_parallel(
  stm_merg,
  #n.samples = 1000,
  cor_method = "pearson",
  remove_zeros = FALSE,
  scale = FALSE,
  log.trans = T,
  log.constant = 1e-6,
  quantile.trim = c(0.01, 0.99),
  use.resampling = F,
  n.cores = parallel::detectCores() - 1
)

pval(corr.pearson.log, use.bonf = T, 0.05)  %>% length()
pval(corr.pearson.log, use.bonf = T, 0.05)
corr_val(corr.pearson.log) #%>% max

pv.sig <- pval(corr.pearson.log, use.bonf = T, 0.05) 
cv.sig <- corr_val(corr.pearson.log)
pairs.sig <- cv.sig[names(pv.sig)]
pairs.sig
pv.sig

# ====================== Correlation: Rate vs. Pair-wise correlation ====================== #
# # Testing whether faster-evolving body parts tend to show stronger correlations 


# ---- cor.test: Pearson ----
pv <- pval(corr.pearson.log, use.bonf = T) 
cv <- corr_val(corr.pearson.log)
cv[names(pv)]


plot(pv, abs(cv))

#Create summed rates for each pair in `cv`
# rate_sums <- sapply(names(cv), function(pair) {
#   parts <- strsplit(pair, " vs ")[[1]]
#   sum(rate_total[parts])
# })

rate_sums <- sapply(names(cv), function(pair) {
  parts <- strsplit(pair, " vs ")[[1]]
  sum(rate_total_norm[parts])
})


plot(rate_sums, abs(cv))
abline(lm(abs(cv) ~ rate_sums), col = "red", lwd = 2)
cor.test(rate_sums, abs(cv), method = "pearson")

# 
# # Filter by pval
# sig_idx <- pv < 0.05
# sig_idx <- pv < 2
# cv_filtered <- cv[sig_idx]
# rate_sums_filtered <- rate_sums[sig_idx]
# 
# #which(rate_sums_filtered<0.0008)
# 
# plot(rate_sums_filtered, abs(cv_filtered))
# abline(lm(abs(cv_filtered) ~ rate_sums_filtered), col = "red", lwd = 2)
# cor.test(rate_sums_filtered, abs(cv_filtered), method = "pearson")
# #cor.test(rate_sums_filtered, abs(cv_filtered), method = "spearman")
# cor.test(rate_sums_filtered, log(abs(cv_filtered)), method = "pearson")
# 
# xy <- trim_xy_quantiles(rate_sums_filtered, abs(cv_filtered), lower = 0.01, upper = 0.99, 
#                   remove_zeros = FALSE, scale = FALSE, log.trans=FALSE, log.constant=0)
# 
# plot(xy$x, xy$y)
# plot(xy$x, log(xy$y))
# cor.test(xy$x, xy$y, method = "pearson")
# cor.test(xy$x, log(xy$y), method = "pearson")

# ---- Zero-Inflated Beta regression: Perhaps we should not report it as it's less data ----
library(gamlss)

pv <- pval(corr.pearson.log, use.bonf = T) 
cv <- corr_val(corr.pearson.log)

# rate_sums <- sapply(names(cv), function(pair) {
#   parts <- strsplit(pair, " vs ")[[1]]
#   sum(rate_total[parts])
# })

rate_sums <- sapply(names(cv), function(pair) {
  parts <- strsplit(pair, " vs ")[[1]]
  sum(rate_total_norm[parts])
})

# Set correlations to 0 if not significant
cv[pv >= 0.05] <- 0  # safer than pv < 0.05 to avoid inversion

# Response: correlation values
cor_vals <- as.numeric(cv)

# Smithson & Verkuilen (2006) transformation to ensure values are inside (0,1)
n <- length(cor_vals)
cor_vals_adj <- (abs(cor_vals) * (n - 1) + 0.5) / n
plot(cor_vals, cor_vals_adj)
cbind(cor_vals, cor_vals_adj)

# Fit Zero-Inflated Beta regression
model_zi_beta <- gamlss(
  formula = cor_vals_adj ~ rate_sums,
  family = BEZI,
  control = gamlss.control(n.cyc = 1000)
)

summary(model_zi_beta)
plot(model_zi_beta)



