# d=1 
# a=1
# b=1e-10 
# t=0.5 
# N=5

hamming_dist <- function(d, a, b, t, N){
  q01 <- (a - a * exp(-1*(a+b) * t) ) / (a+b)
  choose(N, d) * q01^d * (1-q01)^(N-d)
}


# Log-likelihood for full model (a ≠ b)
log_lik_full <- function(params, emp_counts, t, N) {
  a <- params[1]
  b <- params[2]
  if (a <= 0 || b <= 0) return(-Inf)
  
  d_vals <- as.integer(names(emp_counts))
  probs <- sapply(d_vals, function(d) hamming_dist(d, a, b, t, N))
  probs[probs <= 0] <- 1e-12
  sum(emp_counts * log(probs))
}

# Log-likelihood for reduced model (a = b)
log_lik_equal <- function(param, emp_counts, t, N) {
  a <- param[1]
  if (a <= 0) return(-Inf)
  
  d_vals <- as.integer(names(emp_counts))
  probs <- sapply(d_vals, function(d) hamming_dist(d, a, a, t, N))
  probs[probs <= 0] <- 1e-12
  sum(emp_counts * log(probs))
}

compare_models <- function(emp_counts, t, N) {
  # Fit full model (a ≠ b)
  fit_full <- optim(par = c(1, 1e-3),
                    fn = function(p) -log_lik_full(p, emp_counts, t, N),
                    method = "L-BFGS-B", lower = c(1e-5, 1e-10))
  
  # Fit reduced model (a = b)
  fit_equal <- optim(par = c(1),
                     fn = function(p) -log_lik_equal(p, emp_counts, t, N),
                     method = "L-BFGS-B", lower = c(1e-5))
  
  # Extract log-likelihoods
  LL_full <- -fit_full$value
  LL_equal <- -fit_equal$value
  
  # Compute AICs: AIC = -2*logL + 2*k
  AIC_full <- -2 * LL_full + 2 * 2  # two parameters
  AIC_equal <- -2 * LL_equal + 2 * 1  # one parameter
  
  list(
    full = list(params = fit_full$par, logLik = LL_full, AIC = AIC_full),
    equal = list(params = fit_equal$par, logLik = LL_equal, AIC = AIC_equal),
    delta_AIC = AIC_equal - AIC_full
  )
}



log_lik_multi_time <- function(params, timepoints, emp_counts_list, N) {
  a <- params[1]
  b <- params[2]
  if (a <= 0 || b <= 0) return(-Inf)
  
  loglik <- 0
  for (i in seq_along(timepoints)) {
    t <- timepoints[i]
    emp_counts <- emp_counts_list[[i]]
    d_vals <- as.integer(names(emp_counts))
    counts <- as.numeric(emp_counts)
    
    probs <- sapply(d_vals, function(d) hamming_dist(d, a, b, t, N))
    probs[probs <= 0] <- 1e-12
    
    loglik <- loglik + sum(counts * log(probs))
  }
  return(loglik)
}

log_lik_multi_time_equal <- function(param, timepoints, emp_counts_list, N) {
  a <- param[1]
  if (a <= 0) return(-Inf)
  
  loglik <- 0
  for (i in seq_along(timepoints)) {
    t <- timepoints[i]
    emp_counts <- emp_counts_list[[i]]
    d_vals <- as.integer(names(emp_counts))
    counts <- as.numeric(emp_counts)
    
    probs <- sapply(d_vals, function(d) hamming_dist(d, a, a, t, N))
    probs[probs <= 0] <- 1e-12
    
    loglik <- loglik + sum(counts * log(probs))
  }
  return(loglik)
}

compare_multi_time_models <- function(timepoints, emp_counts_list, N) {
  # Fit full model (a ≠ b)
  fit_full <- optim(par = c(1, 0.1),
                    fn = function(p) -log_lik_multi_time(p, timepoints, emp_counts_list, N),
                    method = "L-BFGS-B",
                    lower = c(1e-5, 1e-10))
  
  # Fit reduced model (a == b)
  fit_equal <- optim(par = c(1),
                     fn = function(p) -log_lik_multi_time_equal(p, timepoints, emp_counts_list, N),
                     method = "L-BFGS-B",
                     lower = c(1e-5))
  
  # Extract log-likelihoods
  LL_full <- -fit_full$value
  LL_equal <- -fit_equal$value
  
  # Compute AIC = -2 * logLik + 2 * num_params
  AIC_full <- -2 * LL_full + 2 * 2
  AIC_equal <- -2 * LL_equal + 2 * 1
  
  list(
    full = list(params = fit_full$par, logLik = LL_full, AIC = AIC_full),
    equal = list(params = fit_equal$par, logLik = LL_equal, AIC = AIC_equal),
    delta_AIC = AIC_equal - AIC_full
  )
}

convert_row_to_emp_counts <- function(dh_row) {
  lapply(dh_row, function(d) {
    x <- table(factor(d, levels = 0:max(dh_row, na.rm = TRUE)))
    x[x > 0]
  })
}

#------------

get_state_at_time <- function(map, t_abs) {
  change_times <- unname(map)
  states <- names(map)
  
  if (length(states) == 1) {
    return(states[1])
  }
  
  # Find the first time where t_abs < change time
  idx <- which(t_abs < change_times)[1]
  
  if (is.na(idx)) {
    # If t_abs is beyond all times, return the last state
    return(states[length(states)])
  } else {
    return(states[idx])
  }
}

get_states_all_maps <- function(maps_list, t_abs) {
  sapply(maps_list, get_state_at_time, t_abs = t_abs)
}

get_hamming_distances <- function(maps_list, timepoints) {
  # Helper: get state at a given time
  get_state_at_time <- function(map, t_abs) {
    times <- unname(map)
    states <- names(map)
    if (length(states) == 1) return(states[1])
    idx <- which(t_abs < times)[1]
    if (is.na(idx)) return(states[length(states)])
    return(states[idx])
  }
  
  # Helper: compute Hamming distance
  hamming <- function(s1, s2) {
    sum(strsplit(s1, "")[[1]] != strsplit(s2, "")[[1]])
  }
  
  # Preallocate results: rows = maps, cols = timepoints
  result <- matrix(NA_integer_, nrow = length(maps_list), ncol = length(timepoints))
  colnames(result) <- paste0("t", timepoints)
  
  for (i in seq_along(maps_list)) {
    map <- maps_list[[i]]
    init_state <- names(map)[1]
    for (j in seq_along(timepoints)) {
      current_state <- get_state_at_time(map, timepoints[j])
      result[i, j] <- hamming(init_state, current_state)
    }
  }
  
  as.data.frame(result)
}




#----------------


#Q <- Q_forw
#bit=0
#P_t[bit + 1, ]
simulate_hypercube_ctmc_directional <- function(initial_state, Q_forw, Q_rev, t) {
  if (!all(initial_state %in% c(0, 1))) {
    stop("initial_state must be a vector of 0s and 1s")
  }
  
  library(expm)  # for matrix exponentiation
  
  next_state <- sapply(initial_state, function(bit) {
    Q <- if (bit == 0) Q_forw else Q_rev
    P_t <- expm(Q * t)
    probs <- P_t[bit + 1, ]  # +1 for R indexing
    sample(c(0, 1), size = 1, prob = probs)
  })
  
  return(next_state)
}


simulate_N_hypercube_ctmc_directional <- function(initial_state, Q_forw, Q_rev, t, N) {
  if (!all(initial_state %in% c(0, 1))) {
    stop("initial_state must be a vector of 0s and 1s")
  }
  
  library(expm)  # for matrix exponentiation
  len <- length(initial_state)
  
  replicate(N, {
    sapply(initial_state, function(bit) {
      Q <- if (bit == 0) Q_forw else Q_rev
      P_t <- expm(Q * t)
      probs <- P_t[bit + 1, ]  # +1 because R uses 1-based indexing
      sample(c(0, 1), size = 1, prob = probs)
    })
  }, simplify = "matrix") |> t()
}

hamming_from_initial <- function(init, sim_matrix) {
  if (!is.matrix(sim_matrix)) stop("sim_matrix must be a matrix")
  if (length(init) != ncol(sim_matrix)) stop("init and simulation must have same length")
  
  # Compute Hamming distance row-wise
  apply(sim_matrix, 1, function(row) sum(row != init))
}