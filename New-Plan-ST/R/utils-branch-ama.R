
#' Count state changes per branch in a stochastic map
#'
#' This function takes a stochastic map object (e.g. from 
#' \code{\link[phytools]{make.simmap}}) and returns the number of 
#' state changes that occurred on each branch of the phylogeny.
#'
#' @param simmap A stochastic map object of class \code{"simmap"} 
#'   produced by \code{\link[phytools]{make.simmap}}.
#'
#' @return An integer vector giving the number of state changes per branch, 
#'   in the same order as branches in \code{simmap$maps}.
#'
#' @examples
#' library(phytools)
#' tree <- pbtree(n = 5)
#' states <- sample(c("A","B"), size = 5, replace = TRUE)
#' names(states) <- tree$tip.label
#' simmap <- make.simmap(tree, states, model = "SYM", nsim = 1)
#' count_branch_changes(simmap)
#'
#' @export
count_branch_changes <- function(simmap) {
  changes <- lapply(simmap$maps, function(x) length(x) - 1)
  unlist(changes)
}



#' Summarize branch-specific rates across multiple characters
#'
#' This function takes a list of stochastic maps for multiple characters 
#' (e.g., output from \code{make.simmap}) and computes branch-specific rates 
#' by dividing the number of state changes by branch length. It then sums 
#' these rates across all characters for each stochastic map replicate.
#'
#' @param list.of.chars A list where each element corresponds to a character, 
#'   and contains a list of stochastic map objects (as produced by \code{make.simmap}).
#'
#' @return A list of length equal to the number of stochastic map replicates. 
#'   Each element is a numeric vector giving the summed branch-specific rates 
#'   across all characters for that replicate. 
#'
#' @examples
#'library(phytools)
#'tree <- pbtree(n = 10)
#'states <- sample(1:2, 10, replace = TRUE)
#'names(states) <- tree$tip.label
#'simmap1 <- make.simmap(tree, states, model = "SYM", nsim = 5)
#'simmap2 <- make.simmap(tree, states, model = "SYM", nsim = 5)
#'
#'# Wrap them as a list of characters
#'list.of.chars <- list(char1 = simmap1, char2 = simmap2)
#'
#'# Compute summed rates
#'rates <- ama_bracnhes_multi_char(list.of.chars)
#'print(rates)
#'
#' @export
ama_bracnhes_multi_char <- function(list.of.chars) {
  # Extract branch lengths from the first character, first replicate
  branch.length <- list.of.chars[[1]][[1]]$edge.length
  Nmaps <- length(list.of.chars[[1]])
  
  # Preallocate
  rate.samples <- vector("list", length = length(list.of.chars))
  names(rate.samples) <- names(list.of.chars)
  summed_across_chars <- vector("list", length = Nmaps)
  
  # Convert changes to rates for each character
  for (i in seq_along(list.of.chars)) {
    char <- list.of.chars[[i]]
    rate.samples[[i]] <- lapply(char, function(x) count_branch_changes(x) / branch.length)
  }
  
  # Sum rates across characters per replicate
  for (i in seq_len(Nmaps)) {
    replicate_rates <- lapply(rate.samples, function(x) x[[i]])
    summed_across_chars[[i]] <- Reduce("+", replicate_rates)
  }
  
  return(summed_across_chars)
}


#' Summarize branch-specific rates across multiple characters
#'
#' This function takes a list of stochastic map simulations for multiple 
#' characters and computes summary statistics of branch-specific rates 
#' (e.g., mean, median, standard deviation, or quantiles) across all 
#' replicates and characters.
#'
#' @param list.of.chars A named list of characters, where each element 
#'   contains stochastic map simulations (e.g., from \code{make.simmap}).
#' @param stat Character string specifying the summary statistic to compute. 
#'   Options are \code{"mean"}, \code{"median"}, \code{"sd"}, or \code{"quantile"}.
#' @param probs Numeric vector of probabilities for quantiles (only used if 
#'   \code{stat = "quantile"}). Default is \code{c(0.025, 0.5, 0.975)}.
#'
#' @return A numeric vector or matrix summarizing branch-specific rates:  
#'   \itemize{
#'     \item For \code{"mean"}, \code{"median"}, and \code{"sd"}: a numeric vector 
#'       of length equal to the number of branches.  
#'     \item For \code{"quantile"}: a matrix with rows corresponding to 
#'       quantiles and columns to branches.  
#'   }
#'
#' @examples
#' \dontrun{
#' # Suppose `char_list` is a list of characters with stochastic maps
#' summarize_branch_rates_multi(char_list, stat = "mean")
#' summarize_branch_rates_multi(char_list, stat = "quantile", probs = c(0.1, 0.9))
#' }
#'
#' @seealso \code{\link{summarize_branch_rates}} for single-character summaries.
#'
#' @export
summarize_branch_rates_multi <- function(list.of.chars, 
                                     stat = c("mean", "median", "sd", "quantile"), 
                                     probs = c(0.025, 0.5, 0.975)) {
  stat <- match.arg(stat)
  branch_samples <- ama_bracnhes_multi_char(list.of.chars)
  
  # Convert list of replicates into matrix (rows = replicates, cols = branches)
  mat <- do.call(rbind, branch_samples)
  
  if (stat == "mean") {
    return(colMeans(mat, na.rm = TRUE))
  } else if (stat == "median") {
    return(apply(mat, 2, median, na.rm = TRUE))
  } else if (stat == "sd") {
    return(apply(mat, 2, sd, na.rm = TRUE))
  } else if (stat == "quantile") {
    return(apply(mat, 2, quantile, probs = probs, na.rm = TRUE))
  }
}



#' Summarize branch-specific rates for a single character
#'
#' This function computes summary statistics of branch-specific evolutionary 
#' rates from multiple stochastic map simulations of a single character.
#' Rates are calculated as the number of state changes per branch length, 
#' and then summarized across replicates.
#'
#' @param multi_simmap A list of stochastic map simulations for a single 
#'   character (e.g., from \code{make.simmap}).
#' @param stat Character string specifying the summary statistic to compute. 
#'   Options are \code{"mean"}, \code{"median"}, \code{"sd"}, or \code{"quantile"}.
#' @param probs Numeric vector of probabilities for quantiles (only used if 
#'   \code{stat = "quantile"}). Default is \code{c(0.025, 0.5, 0.975)}.
#'
#' @return A numeric vector or matrix summarizing branch-specific rates:  
#'   \itemize{
#'     \item For \code{"mean"}, \code{"median"}, and \code{"sd"}: a numeric vector 
#'       of length equal to the number of branches.  
#'     \item For \code{"quantile"}: a matrix with rows corresponding to 
#'       quantiles and columns to branches.  
#'   }
#'
#' @examples
#' \dontrun{
#' # Suppose `char_simmap` is a list of stochastic map replicates for one character
#' summarize_branch_rates(char_simmap, stat = "mean")
#' summarize_branch_rates(char_simmap, stat = "quantile", probs = c(0.1, 0.9))
#' }
#'
#' @seealso \code{\link{summarize_branch_rates_multi}} for summarizing across multiple characters.
#'
#' @export
summarize_branch_rates <- function(multi_simmap, 
                                   stat = c("mean", "median", "sd", "quantile"), 
                                   probs = c(0.025, 0.5, 0.975)) {
  stat <- match.arg(stat)
  branch.length <- multi_simmap[[1]]$edge.length
  branch_samples <- lapply(multi_simmap, function(x) count_branch_changes(x) / branch.length)
  
  # Convert list of replicates into matrix (rows = replicates, cols = branches)
  mat <- do.call(rbind, branch_samples)
  
  if (stat == "mean") {
    return(colMeans(mat, na.rm = TRUE))
  } else if (stat == "median") {
    return(apply(mat, 2, median, na.rm = TRUE))
  } else if (stat == "sd") {
    return(apply(mat, 2, sd, na.rm = TRUE))
  } else if (stat == "quantile") {
    return(apply(mat, 2, quantile, probs = probs, na.rm = TRUE))
  }
}