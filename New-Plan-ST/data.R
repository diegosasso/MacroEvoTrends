# Load libraries.
library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)
source("R/utils-ST.R")

# ====================== Data ====================== #
# Data:
# - 15 body regions (first order)
# - 5 second level BD 
# - 1 sclerites dataset, and 1 muscles  dataset
# - 1 larval dataset 


# ---- Number of states per BR ----
load('data/char_num.RDA')
anat_ent_state
# log of states per character per BR
log_states <- lapply(anat_ent_state, function(x) sum(log(get_len(x))) ) %>% unlist
tot_states <- lapply(anat_ent_state, function(x) sum(get_len(x)) ) %>% unlist

# ---- 15 body regions (first order) ----
# Read in stm_merge
stm_merg <- readRDS("data/stm_merg.RDS")
names(stm_merg)

# ---- 5 second level BD  ----
stm_merg2 <- c()
stm_merg2$'head' <- do.call(c, stm_merg[1:2])
stm_merg2$'mesosoma' <- do.call(c, stm_merg[3:8])
stm_merg2$'metasoma_gen' <- do.call(c, stm_merg[9:10])
stm_merg2$'legs' <- do.call(c, stm_merg[11:13])
stm_merg2$'wings' <- do.call(c, stm_merg[14:15])
stm_merg2

# saveRDS(stm_merg2, "data/stm_merg2.RDS")

# ---- 5 second level BD  ----
stm_merg2 <-  readRDS("data/stm_merg2.RDS")
stm_merg2
names(stm_merg2)
rates_list2 <- lapply(stm_merg2, get_branch_rate_across_smaps)
rates2 <- do.call(cbind, rates_list2)
colMeans(rates2 == 0)
names(stm_merg2)

# ====================== Note:  Branch specific rate ====================== #
smap <- stm_merg$cranium[[1]]
# rate for branch 151
rate_151 <- length(smap$maps[[151]]-1)/smap$edge.length[[151]]
rate_151