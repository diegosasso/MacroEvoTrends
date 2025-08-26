library(expm)
library(ape)
library(phytools)
library(corHMM)
library(rphenoscate)
library(ontophylo)
library(dplyr)
source('R/lik-hamming.R')


tree <- readRDS("tree_test.RDS")
# plot.phylo(tree, cex = 0.5)
# nodelabels(frame = "none", col = "blue", cex = 0.5)
# edgelabels(frame = "none", col = "blue", cex = 0.5)


tree.paint <- paintBranches(tree, edge=167, state=1, anc.state="0")
plot(tree.paint)



Q_base <- initQ(c(0, 1), c(0.01,0.01))
Q_base
Q_jump <- initQ(c(0, 1), c(10, 1e-10))
Q_jump
Q_list <- list('0'=Q_base, '1'=Q_jump)

char <- sim.multiMk(tree.paint, Q_list, anc='0', nsim=2)
str(char)
#-------
Q_jump
expm(Q_jump*0.1)

#-----

Q_est <- initQ(c(0, 1), c(1,1), diag.as = NA)
Q_est.list <- list('0'=Q_est, '1'=Q_est)
fitmultiMk(tree.paint, char, model="ARD")
fitmultiMk(tree.paint, char, model=Q_est)

model = 'ARD'

tree <- tree.paint
x <- char
x <- char[,1]
fitmultiMk_mine(tree.paint, char)
fitmultiMk(tree.paint, char, model="ER")

fitmultiMk_mine <- function (tree, x) 
{
  if (!inherits(tree, "simmap")) {
    stop("tree should be an object of class \"simmap\". Use fitMk.\n")
  }
  else {
    regimes <- mapped.states(tree)
    nregimes <- length(regimes)
  }
  if (hasArg(q.init)) 
    q.init <- list(...)$q.init
  else q.init <- length(unique(x))/sum(tree$edge.length)
  if (hasArg(rand_start)) 
    rand_start <- list(...)$rand_start
  else rand_start <- FALSE
  if (hasArg(opt.method)) 
    opt.method <- list(...)$opt.method
  else opt.method <- "nlminb"
  if (hasArg(min.q)) 
    min.q <- list(...)$min.q
  else min.q <- 1e-12
  # if (is.matrix(x)) {
  #   x <- x[tree$tip.label, ]
  #   m <- ncol(x)
  #   states <- colnames(x)
  # }
  # else {
  #   x <- to.matrix(x, sort(unique(x)))
  #   x <- x[tree$tip.label, ]
  #   m <- ncol(x)
  #   states <- colnames(x)
  # }
  #
  #---- make characters
  X <- vector(mode = 'list', length = ncol(x))
  i=1
  for (i in 1:ncol(x)){
    X.i <-  to.matrix(x[,i], sort(unique(x[,i])))
    X.i <- X.i[tree$tip.label, ]
    X[[i]] <- X.i
  }
  m <- 2
  states <- c("0", "1")
  #--------------
  if (hasArg(pi)) 
    pi <- list(...)$pi
  else pi <- "equal"
  if (pi[1] == "equal") 
    pi <- setNames(rep(1/m, m), states)
  else pi <- pi/sum(pi)
  #
  # if (is.character(model)) {
  #   rate <- matrix(NA, m, m)
  #   if (model == "ER") {
  #     k <- rate[] <- 1
  #     diag(rate) <- NA
  #   }
  #   else if (model == "ARD") {
  #     k <- m * (m - 1)
  #     rate[col(rate) != row(rate)] <- 1:k
  #   }
  #   else if (model == "SYM") {
  #     k <- m * (m - 1)/2
  #     ii <- col(rate) < row(rate)
  #     rate[ii] <- 1:k
  #     rate <- t(rate)
  #     rate[ii] <- 1:k
  #   }
  # }
  # else {
  #   if (ncol(model) != nrow(model)) 
  #     stop("model is not a square matrix")
  #   if (ncol(model) != ncol(x)) 
  #     stop("model does not have the right number of columns")
  #   rate <- model
  #   k <- max(rate)
  # }
  #rate <- matrix(NA, m, m)
  Q <- replicate(nregimes, matrix(0, m, m), simplify = FALSE)
  names(Q) <- regimes
  #index.matrix <- rate
  tmp <- cbind(1:m, 1:m)
  #rate[tmp] <- 0
  #rate[rate == 0] <- k + 1
  pw <- reorder(map.to.singleton(tree), "pruningwise")
  N <- Ntip(pw)
  M <- pw$Nnode
  
  liks <- rbind(x, matrix(0, M, m, dimnames = list(1:M + N, 
                                                   states)))
  #------
  # pp <- rep(0.08, 4)
  # output.liks = FALSE
  #-----
  lik <- function(pp, output.liks = FALSE, pi) {
    if (any(is.nan(pp)) || any(is.infinite(pp))) 
      return(1e+50)
    comp <- vector(length = N + M, mode = "numeric")
    # for (i in 1:nregimes) {
    #   Q[[i]][] <- c(pp[1:k + (i - 1) * k], 0)[rate]
    #   diag(Q[[i]]) <- -rowSums(Q[[i]])
    # }
    #----------
    Q[[1]] <- rphenoscate::initQ(c(0, 1), c(pp[1],pp[1]))
    Q[[2]] <- rphenoscate::initQ(c(0, 1), c(pp[2],pp[3]))
    #----------
    parents <- unique(pw$edge[, 1])
    root <- min(parents)
    i=1
    for (i in 1:length(parents)) {
      anc <- parents[i]
      ii <- which(pw$edge[, 1] == parents[i])
      desc <- pw$edge[ii, 2]
      el <- pw$edge.length[ii]
      v <- vector(length = length(desc), mode = "list")
      reg <- names(pw$edge.length)[ii]
      for (j in 1:length(v)) v[[j]] <- phytools:::EXPM(Q[[reg[j]]] * 
                                              el[j]) %*% liks[desc[j], ]
      vv <- if (anc == root) 
        Reduce("*", v)[, 1] * pi
      else Reduce("*", v)[, 1]
      comp[anc] <- sum(vv)
      liks[anc, ] <- vv/comp[anc]
    }
    logL <- -sum(log(comp[1:M + N]))
    return(if (is.na(logL)) Inf else logL)
  }
  #------
  # pp <- rep(0.08, 4)
  # output.liks = FALSE
  #-----
  lik_one <- function(pp, output.liks = FALSE, pi) {
    if (any(is.nan(pp)) || any(is.infinite(pp))) 
      return(1e+50)
    comp <- vector(length = N + M, mode = "numeric")
    # for (i in 1:nregimes) {
    #   Q[[i]][] <- c(pp[1:k + (i - 1) * k], 0)[rate]
    #   diag(Q[[i]]) <- -rowSums(Q[[i]])
    # }
    #----------
    Q[[1]] <- rphenoscate::initQ(c(0, 1), c(pp[1],pp[1]))
    Q[[2]] <- rphenoscate::initQ(c(0, 1), c(pp[2],pp[3]))
    #----------
    parents <- unique(pw$edge[, 1])
    root <- min(parents)
    i=1
    for (i in 1:length(parents)) {
      anc <- parents[i]
      ii <- which(pw$edge[, 1] == parents[i])
      desc <- pw$edge[ii, 2]
      el <- pw$edge.length[ii]
      v <- vector(length = length(desc), mode = "list")
      reg <- names(pw$edge.length)[ii]
      for (j in 1:length(v)) v[[j]] <- phytools:::EXPM(Q[[reg[j]]] * 
                                                         el[j]) %*% liks[desc[j], ]
      vv <- if (anc == root) 
        Reduce("*", v)[, 1] * pi
      else Reduce("*", v)[, 1]
      comp[anc] <- sum(vv)
      liks[anc, ] <- vv/comp[anc]
    }
    logL <- -sum(log(comp[1:M + N]))
    return(if (is.na(logL)) Inf else logL)
  }
  ####
  # if (length(q.init) != (k * nregimes)) 
  #   q.init <- rep(q.init[1], k * nregimes)
  if (length(q.init) == 1) 
    q.init <- rep(q.init[1], 3)
  if (rand_start) 
    q.init <- q.init * rexp(length(q.init), 1)
  if (opt.method == "optim") 
    fit <- optim(q.init, function(p) lik(p, pi = pi), method = "L-BFGS-B", 
                 lower = rep(min.q, 3))
  else fit <- nlminb(q.init, function(p) lik(p, pi = pi), 
                     lower = rep(0, 3), upper = rep(1e+50, 3))
  Q[[1]] <- rphenoscate::initQ(c(0, 1), c(fit$par[1],fit$par[1]))
  Q[[2]] <- rphenoscate::initQ(c(0, 1), c(fit$par[2],fit$par[3]))
  obj <- list(logLik = if (opt.method == "optim") -fit$value else -fit$objective, 
              rates = Q, #index.matrix = index.matrix, 
              states = states, 
              regimes = regimes, pi = pi, method = opt.method)
  lik.f <- function(q) -lik(q, output.liks = FALSE, pi = pi)
  obj$lik <- lik.f
  #class(obj) <- "fitmultiMk"
  return(obj)
}

