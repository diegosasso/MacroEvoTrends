library(expm)
source('R/lik-hamming.R')

Q_forw <- initQ(c(0, 1), c(1, 1e-10))
Q_rev <- initQ(c(0, 1), c(1e-10, 1))


# Simulate
init <- c(1, 0, 1, 0, 1)
init2 <- c(0, 0, 0, 0, 0)
t <- 0.5

#simulate_hypercube_ctmc_directional(init, Q_forw, Q_rev, t)
sim_matrix <- simulate_N_hypercube_ctmc_directional(init, Q_forw, Q_rev, t, N=10000)
sim_matrix2 <- simulate_N_hypercube_ctmc_directional(init2, Q_forw, Q_rev, t, N=10000)

Dh <- hamming_from_initial(init, sim_matrix)
Dh2 <- hamming_from_initial(init2, sim_matrix2)
hist(Dh, probability = TRUE)
hist(Dh2, probability = T, col = rgb(0, 0, 1, 0.3), add = TRUE)

#------- Plot Emp vs. Analit   probabilities

# Your empirical distribution
emp_freq <- table(Dh)
emp_prob <- prop.table(emp_freq)

# Analytical probability
dist <- 0:5
Pr <- sapply(dist, function(X) hamming_dist(d=X, a=1, b=1e-10, t=0.5, N=5))

# Plot   probabilities
bp <- barplot(emp_prob,
              names.arg = dist,
              ylim = c(0, max(c(emp_prob, Pr)) + 0.05),
              col = "skyblue",
              xlab = "Hamming Distance",
              ylab = "Probability",
              main = "Empirical vs Analytical Hamming Distribution")
points(x = bp, y = Pr, type = "b", pch = 16, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "Analytical"),
       fill = c("skyblue", NA), border = NA,
       lty = c(NA, 1), pch = c(NA, 16), col = c("skyblue", "red"))




#------ Likelihood

Q_forw <- initQ(c(0, 1), c(1, 1e-10))
Q_rev <- initQ(c(0, 1), c(1e-10, 1))

# Parameters
t <- 2
N <- 30

init <- rep(0, N)
sim_matrix <- simulate_N_hypercube_ctmc_directional(init, Q_forw, Q_rev, t, N=1000)
Dh <- hamming_from_initial(init, sim_matrix)
hist(Dh)
emp_counts <- table(Dh)



# Optimize (maximize log-likelihood → minimize negative log-likelihood)
optim_result <- optim(par = c(1, 1e-3), 
                      fn = function(p) -log_lik_full(p, emp_counts, t, N),
                      method = "L-BFGS-B", lower = c(1e-5, 1e-10))

optim_result$par  # estimated a and b


compare_models(emp_counts, t, N)





#------ Multi-timepoint likelihood fitting

library(expm)
library(ape)
library(phytools)
library(corHMM)
library(rphenoscate)
library(ontophylo)
library(dplyr)
source('R/lik-hamming.R')

# Simulate empirical data at multiple time points

timepoints <- c(0.1, 1, 3)
N_bits <- 30
N_samples <- 1

Q_forw <- initQ(c(0, 1), c(.5, 1e-10))
Q_rev <- initQ(c(0, 1), c(1e-10, .5))
init <- rep(0, N_bits)

emp_counts_list <- lapply(timepoints, function(t) {
  sim_matrix <- simulate_N_hypercube_ctmc_directional(init, Q_forw, Q_rev, t, N = N_samples)
  Dh <- hamming_from_initial(init, sim_matrix)
  table(Dh)
})

emp_counts_list

neg_loglik <- function(params) {
  -log_lik_multi_time(params, timepoints, emp_counts_list, N = N_bits)
}

optim_result <- optim(par = c(1, 0.1), fn = neg_loglik,
                      method = "L-BFGS-B", lower = c(1e-5, 1e-10))

optim_result$par  # Estimated a and b

# ----- Compare Models

result <- compare_multi_time_models(timepoints, emp_counts_list, N_bits)
result
result$full$params     # estimated a, b
result$equal$params    # estimated a (where a == b)
result$delta_AIC       # < 0 → full model better


#--------- Real data

library(expm)
library(ape)
library(phytools)
library(corHMM)
library(rphenoscate)
library(ontophylo)
library(dplyr)
source('R/lik-hamming.R')

stm_merg_jump <- readRDS('RDS/stm_merg_toy-ST.RDS')
# stm_merg_jump <- readRDS('RDS/stm_merg_jump-ST.RDS')
# stm_merg_sym <- readRDS('RDS/stm_merg_sym-ST.RDS')
#stm_merg_toy <- readRDS('RDS/stm_merg_toy-ST.RDS')
tree <- readRDS("tree_test.RDS")
plot.phylo(tree, cex = 0.5)
edgelabels(frame = "none", col = "blue", cex = 0.5)

br_maps <- lapply(stm_merg_jump, function(x) x$maps[[188]] )
lapply(br_maps, length) %>% unlist
br_maps_cumsum <- lapply(br_maps, cumsum)

# Extract Hamming distances
branch_length <- stm_merg_jump[[1]]$edge.length[188]
#timepoints <- c(0.01, 0.02, 0.03, 0.05, 0.075)
timepoints <- seq(0.001, branch_length, length.out = 10)
Dh_matrix <- get_hamming_distances(br_maps_cumsum, timepoints)

head(br_maps_cumsum)
head(Dh_matrix)
emp_counts_list <- convert_row_to_emp_counts(Dh_matrix[5,])

N_bits=50
result <- compare_multi_time_models(timepoints, emp_counts_list, N_bits)
result

N_bits=50
branch_length <- stm_merg_jump[[1]]$edge.length[188]
timepoints <- seq(0.001, branch_length, length.out = 500)
Dh_matrix <- get_hamming_distances(br_maps_cumsum, timepoints)
dAIC <- c()
i=1
for (i in 1:100){
  emp_counts_list <- convert_row_to_emp_counts(Dh_matrix[i,])
  result <- compare_multi_time_models(timepoints, emp_counts_list, N_bits)
  dAIC <- c(dAIC, result$delta_AIC)
}
dAIC
hist(dAIC)
sum(dAIC > 2)
