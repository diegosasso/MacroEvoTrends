#---------
# cannot use ?
# check chars to be non invariant
# in PARAMO, do not use root yang
#----------

library(expm)
library(ape)
library(phytools)
library(corHMM)
library(rphenoscate)
library(ontophylo)
library(dplyr)
source('R-hmm/utils-hmm.R')
source('R-hmm/dev.rayDisc-multi.R')
# 
# tree <- readRDS("tree_test.RDS")
# char <- readRDS('RDS/char.rds')
# 
# phy <- tree
# # data <- cbind(rownames(char), char[,1]+1)
# data <- cbind(rownames(char), char)
# str(data)
# Nchar = 2
# 
# phy
# data
# #ntraits = 1
# #charnum = 1
# #rate.mat = NULL
# rate.mat = initQ(c(0, 1), c(1, 1), diag.as = NA)
# #model = c("ER", "SYM", "ARD"),
# node.states = c("none")
# state.recon = c("subsequently")
# lewis.asc.bias = FALSE
# #p = NULL
# p=c(0.1, 0.1)
# root.p = "flat"
# ip = NULL
# lb = 1e-09
# ub = 100
# verbose = TRUE
# diagn = FALSE
# hmm.map <-  c("0&1", "2&3")




rayDISC_multi <- function (phy, data, Nchar, rate.mat = NULL, hmm.map, add.dummy=FALSE,
          node.states = "none", state.recon = c("subsequently"), 
          lewis.asc.bias = FALSE, p = NULL, root.p = "flat", ip = NULL, 
          lb = 1e-09, ub = 100, verbose = TRUE, diagn = FALSE) 
{
  if (is.null(node.states)) {
    obj <- NULL
    obj$loglik <- NULL
    obj$diagnostic <- paste("No model for ancestral states selected.  Please pass one of the following to rayDISC command for parameter 'node.states': joint, marginal, or scaled.")
    return(obj)
  }else {
    valid.models <- c("joint", "marginal", "scaled", "none")
    if (!any(valid.models == node.states)) {
      obj <- NULL
      obj$loglik <- NULL
      obj$diagnostic <- paste("'", node.states, "' is not valid for ancestral state reconstruction method.  Please pass one of the following to rayDISC command for parameter 'node.states': joint, marginal, or scaled.", 
                              sep = "")
      return(obj)
    }
    if (length(node.states) > 1) {
      node.states <- "marginal"
      cat("No model selected for 'node.states'. Will perform marginal ancestral state estimation.\n")
    }
  }
  #
  if (!state.recon == "subsequently" & node.states == "marginal" | 
      node.states == "scaled") {
    stop("Simultaneous estimation of rates and states using either marginal or scaled probabilities not yet implemented.", 
         call. = FALSE)
  }
  if (!state.recon == "subsequently") {
    if (!is.null(phy$node.label)) {
      if (!is.na(phy$node.label[Ntip(phy) + 1])) {
        root.p <- NULL
      }
    }
  }
  if (!state.recon == "estimate" & !state.recon == "given" & 
      !state.recon == "subsequently") {
    stop("Check that you have a supported state.recon analysis. Options are subsequently, estimate, or given.", 
         call. = FALSE)
  }
  if (state.recon == "subsequently") {
    phy$node.label <- NULL
  }else {
    if (state.recon == "given") {
      if (is.null(phy$node.label)) {
        stop("You specified you wanted to estimate rates on a given character history, but the tree does not contain node labels.", 
             call. = FALSE)
      }else {
        if (any(is.na(phy$node.label))) {
          cat("Model will assume you want to estimate rates and states, but include state constraints on some but not all nodes.\n")
        }else {
          cat("Model will assume you want to estimate rates, but include state constraints all nodes.\n")
        }
      }
    }else {
      if (is.null(phy$node.label)) {
        cat("Model will assume you want to estimate rates and states simultaneously.\n")
      }else {
        cat("Model will assume you want to estimate rates and states, but include state constraints on some but not all nodes.\n")
        state.recon = "given"
      }
    }
  }
  if (!is.null(root.p)) {
    if (!is.character(root.p)) {
      root.p <- root.p/sum(root.p)
    }
  }
  phy$edge.length[phy$edge.length == 0] = 1e-05
  matching <- corHMM:::match.tree.data(phy, data)
  data <- matching$data
  phy <- matching$phy
  #--
  # data[(data[, charnum + 1] == "?"), charnum + 1] <- NA
  # if (nlevels(as.factor(data[, charnum + 1])) <= 1) {
  #   obj <- NULL
  #   obj$loglik <- NULL
  #   obj$diagnostic <- paste("Character ", charnum, " is invariant. Analysis stopped.", 
  #                           sep = "")
  #   return(obj)
  # }else {
  #   lvls <- as.factor(data[, charnum + 1])
  #   if (nlevels(as.factor(data[, charnum + 1])) == 2 && 
  #       any(lvls %in% c("?", "NA"))) {
  #     obj <- NULL
  #     obj$loglik <- NULL
  #     obj$diagnostic <- paste("Character ", charnum, " is invariant. Analysis stopped.", 
  #                             sep = "")
  #     return(obj)
  #   }
  # }
  # charnum=1
  # add.dummy is when the char does not contain all states
  char_env <-prepare_vectorized_chars(data, phy, Nchar)
  char_env_summary(char_env)
  # make invar chars
  char_invar_env <-make_invariant_env(data, phy, hmm.map)
  #ls(char_env)
  #char_env$char2
  #
  #-- this is needed for model.set.final <- corHMM:::rate.cat.set.rayDISC
  data.sort <- data.frame(data[, 1 + 1], data[, 1 + 1], row.names = data[, 1])
  data.sort <- data.sort[phy$tip.label, ]
  #-------
  # this is a quick fix if character does not contian all states
  if (add.dummy){
    dummy <- c(paste0(hmm.map, collapse = '&'), paste0(hmm.map, collapse = '&'))
    data.sort.mod <- data.sort
    data.sort.mod <-rbind(dummy, data.sort.mod)
    model.set.final_D <- corHMM:::rate.cat.set.rayDISC(phy = phy, data = data.sort.mod, model = 'ER', charnum = 1)
    model.set.final_D$liks <- model.set.final_D$liks[-1,]
    char_env[['char1']] <- model.set.final_D$liks
  }
  # data.sort <- data.frame(data[, charnum + 1], data[, charnum + 1], row.names = data[, 1])
  # data.sort <- data.sort[phy$tip.label, ]
  # data.rayDISC <- data.frame(sp = rownames(data.sort), d = data.sort[,1])
  # counts <- table(data.sort[, 1])
  # levels <- levels(as.factor(data.sort[, 1]))
  # cols <- as.factor(data.sort[, 1])
  # if (verbose == TRUE) {
  #   cat("State distribution in data:\n")
  #   cat("States:", levels, "\n", sep = "\t")
  #   cat("Counts:", counts, "\n", sep = "\t")
  # }
  # k <- 1
  # factored <- corHMM:::factorData(data.sort, charnum = charnum)
  #-------
  # nl <- ncol(factored)
  # state.names <- colnames(factored)
  nl <- ncol(rate.mat)
  state.names <- colnames(rate.mat)
  bound.hit <- FALSE
  if (ub < 0) {
    ub <- log(100)
  }else {
    ub <- log(ub)
  }
  if (lb <= 0) {
    lb <- -21
  }else {
    lb <- -21
  }
  if (ub < lb) {
    ub <- log(100)
    lb <- -21
  }
  obj <- NULL
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  #model = model
  root.p = root.p
  ip = ip
  # model.set.final <- corHMM:::rate.cat.set.rayDISC(phy = phy, data = data.sort, 
  #                                         model = model, charnum = charnum)
  # model.set.final <-corHMM:::rate.cat.set.rayDISC(phy = phy, data = data.sort, model = 'ARD', charnum = 1)
  # model.set.final$rate
  model.set.final <- corHMM:::rate.cat.set.rayDISC(phy = phy, data = data.sort, model = 'ER', charnum = 1)
  if (add.dummy){
    model.set.final$Q <- model.set.final_D$Q
    model.set.final$rate <-  model.set.final_D$rate
  }
  if (!is.null(rate.mat)) {
    rate <- rate.mat
    model.set.final$np <- max(rate, na.rm = TRUE)
    rate[is.na(rate)] = max(rate, na.rm = TRUE) + 1
    model.set.final$rate <- rate
    model.set.final$index.matrix <- rate.mat
  }
  lower = rep(lb, model.set.final$np)
  upper = rep(ub, model.set.final$np)
  opts <- list(algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000", 
               ftol_rel = .Machine$double.eps^0.5)
  if (!is.null(p)) {
    if (verbose == TRUE) {
      cat("Calculating likelihood from a set of fixed parameters", 
          "\n")
    }
    out <- NULL
    out$solution <- p
    phy <- reorder(phy, "pruningwise")
    if (state.recon == "subsequently") {
      # out$objective <- dev.raydisc(log(out$solution), 
      #                              phy = phy, liks = model.set.final$liks, Q = model.set.final$Q, 
      #                              rate = model.set.final$rate, root.p = root.p, 
      #                              lewis.asc.bias = lewis.asc.bias)
      out$objective <- dev.raydisc_multi(
                                  log(out$solution), 
                                  phy = phy, char_env = char_env, char_invar_env=char_invar_env, Q = model.set.final$Q, 
                                  rate = model.set.final$rate, root.p = root.p, 
                                  lewis.asc.bias = lewis.asc.bias)
      loglik <- -out$objective
    }
    est.pars <- out$solution
  }
  else {
    if (is.null(ip)) {
      if (verbose == TRUE) {
        cat("Initializing...", "\n")
      }
      model.set.init <- corHMM:::rate.cat.set.rayDISC(phy = phy, data = data.sort, model = "ER", charnum = 1)
      opts <- list(algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000", ftol_rel = .Machine$double.eps^0.5)
      
      par.score_vec <- c()
      for (charnum in 1:Nchar){
        data.sort <- data.frame(data[, charnum + 1], data[, charnum + 1], row.names = data[, 1])
        data.sort <- data.sort[phy$tip.label, ]
        taxa.missing.data.drop <- which(is.na(data.sort[, 1]))
        if (length(taxa.missing.data.drop) != 0) {
          tip.labs <- names(taxa.missing.data.drop)
          dat <- as.matrix(data.sort)
          dat.red <- dat[-taxa.missing.data.drop, ]
          phy.red <- drop.tip(phy, taxa.missing.data.drop)
          dat.red <- phangorn::phyDat(dat.red, type = "USER", levels = sort(unique(c(dat))))
          par.score <- phangorn::parsimony(phy.red, dat.red, method = "fitch")/2
          par.score_vec <- c(par.score_vec, par.score)
        }else {
          dat <- as.matrix(data.sort)
          dat <- phangorn::phyDat(dat, type = "USER", levels = sort(unique(c(dat))))
          phy.tmp <- ape::multi2di(phy)
          par.score <- phangorn::parsimony(phy.tmp, dat, method = "fitch")/2
          par.score_vec <- c(par.score_vec, par.score)
        }
      }

      tl <- sum(phy$edge.length)
      # mean.change = par.score/tl
      mean.change = sum(par.score_vec)/tl
      if (mean.change == 0) {
        ip = 0.01
      }else {
        ip <- rexp(1, 1/mean.change)
      }
      if (log(ip) < lb || log(ip) > ub) {
        ip <- exp(lb)
      }
      lower.init = rep(lb, model.set.init$np)
      upper.init = rep(ub, model.set.init$np)
      phy <- reorder(phy, "pruningwise")
      # init = nloptr(x0 = rep(log(ip), length.out = model.set.init$np), 
      #               eval_f = dev.raydisc, lb = lower.init, ub = upper.init, 
      #               opts = opts, phy = phy, liks = model.set.init$liks, 
      #               Q = model.set.init$Q, rate = model.set.init$rate, 
      #               root.p = root.p, lewis.asc.bias = lewis.asc.bias)
      init = nloptr(x0 = rep(log(ip), length.out = model.set.init$np), 
                    eval_f = dev.raydisc_multi, lb = lower.init, ub = upper.init, 
                    opts = opts, phy = phy, char_env=char_env, char_invar_env=char_invar_env,
                    Q = model.set.init$Q, rate = model.set.init$rate, 
                    root.p = root.p, lewis.asc.bias = lewis.asc.bias)
      if (verbose == TRUE) {
        cat("Finished. Beginning thorough search...", 
            "\n")
      }
      lower = rep(lb, model.set.final$np)
      upper = rep(ub, model.set.final$np)
      if (state.recon == "subsequently") {
        out <- nloptr(x0 = rep(init$solution, length.out = model.set.final$np), 
                      eval_f = dev.raydisc_multi, lb = lower, ub = upper, 
                      opts = opts, phy = phy, char_env=char_env, char_invar_env=char_invar_env,
                      Q = model.set.final$Q, rate = model.set.final$rate, 
                      root.p = root.p, lewis.asc.bias = lewis.asc.bias)
      }
      # else {
      #   out <- nloptr(x0 = rep(init$solution, length.out = model.set.final$np), 
      #                 eval_f = dev.raydisc.rates.and.states, lb = lower, 
      #                 ub = upper, opts = opts, phy = phy, data = data, 
      #                 hrm = FALSE, rate.cat = NULL, rate.mat = rate.mat, 
      #                 ntraits = ntraits, method = node.states, model = model, 
      #                 charnum = charnum, root.p = root.p, lewis.asc.bias = lewis.asc.bias, 
      #                 get.likelihood = TRUE)
      # }
      loglik <- -out$objective
      est.pars <- exp(out$solution)
    #--------- ip is givem
    } else {
      phy <- reorder(phy, "pruningwise")
      if (verbose == TRUE) {
        cat("Beginning subplex optimization routine -- Starting value(s):",
            ip, "\n")
      }
      opts <- list(algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000",
                   ftol_rel = .Machine$double.eps^0.5)
      if (state.recon == "subsequently") {
        if (!length(ip) == model.set.final$np)
          stop(" Length of starting state vector does not match model parameters. ")
        out <- nloptr(x0 = log(ip), eval_f = dev.raydisc_multi,
                      lb = lower, ub = upper, opts = opts, phy = phy,
                      char_env=char_env, char_invar_env=char_invar_env, Q = model.set.final$Q,
                      rate = model.set.final$rate, root.p = root.p,
                      lewis.asc.bias = lewis.asc.bias)
      }
      # else {
      #   if (!length(ip) == model.set.final$np)
      #     stop(" Length of starting state vector does not match model parameters. ")
      #   out <- nloptr(x0 = log(ip), eval_f = dev.raydisc.rates.and.states,
      #                 lb = lower, ub = upper, opts = opts, phy = phy,
      #                 data = data, hrm = FALSE, rate.cat = NULL,
      #                 rate.mat = rate.mat, ntraits = ntraits, method = node.states,
      #                 model = model, charnum = charnum, root.p = root.p,
      #                 lewis.asc.bias = lewis.asc.bias, get.likelihood = TRUE)
      # }
      loglik <- -out$objective
      est.pars <- exp(out$solution)
    }
  }
  if (verbose == TRUE) {
    cat("Finished. Inferring ancestral states using", node.states, 
        "reconstruction.", "\n")
  }
  TIPS <- 1:nb.tip
  if (node.states == "none") {
    lik.anc <- NULL
    lik.anc$lik.tip.states <- "You turned this feature off. Try plugging into ancRECON function directly."
    lik.anc$lik.anc.states <- "You turned this feature off. Try plugging into ancRECON function directly."
    tip.states <- lik.anc$lik.tip.states
  }
  # else {
  #   data.rayDISC[is.na(data.rayDISC)] <- "?"
  #   if (node.states == "marginal" || node.states == "scaled") {
  #     lik.anc <- ancRECON(phy, data.rayDISC, est.pars, 
  #                         rate.cat = NULL, rate.mat = rate.mat, ntraits = ntraits, 
  #                         method = node.states, model = model, root.p = root.p)
  #     pr <- apply(lik.anc$lik.anc.states, 1, which.max)
  #     phy$node.label <- pr
  #     tip.states <- lik.anc$lik.tip.states
  #   }
  #   if (!state.recon == "given") {
  #     if (node.states == "joint") {
  #       lik.anc <- ancRECON(phy, data.rayDISC, est.pars, 
  #                           rate.cat = NULL, rate.mat = rate.mat, ntraits = ntraits, 
  #                           method = node.states, model = model, root.p = root.p)
  #       phy$node.label <- lik.anc$lik.anc.states
  #       tip.states <- lik.anc$lik.tip.states
  #     }
  #   }
  #   else {
  #     if (any(is.na(phy$node.label))) {
  #       lik.anc <- ancRECON(phy, data.rayDISC, est.pars, 
  #                           rate.cat = NULL, rate.mat = rate.mat, ntraits = ntraits, 
  #                           method = node.states, model = model, root.p = root.p)
  #       phy$node.label <- lik.anc$lik.anc.states
  #       tip.states <- lik.anc$lik.tip.states
  #     }
  #     else {
  #       lik.anc <- NULL
  #       lik.anc$lik.anc.states <- phy$node.label
  #       lik.anc$lik.tip.states <- data.sort[, 1]
  #       tip.states <- lik.anc$lik.tip.states
  #     }
  #   }
  # }
  if (diagn == TRUE) {
    if (verbose == TRUE) {
      cat("Finished. Performing diagnostic tests.", "\n")
    }
    h <- hessian(func = dev.raydisc, x = log(est.pars), 
                 phy = phy, liks = model.set.final$liks, Q = model.set.final$Q, 
                 rate = model.set.final$rate, root.p = root.p, lewis.asc.bias = lewis.asc.bias)
    solution <- matrix(est.pars[model.set.final$index.matrix], 
                       dim(model.set.final$index.matrix))
    solution.se <- matrix(sqrt(diag(pseudoinverse(h)))[model.set.final$index.matrix], 
                          dim(model.set.final$index.matrix))
    hess.eig <- eigen(h, symmetric = TRUE)
    eigval <- signif(hess.eig$values, 2)
    eigvect <- round(hess.eig$vectors, 2)
  }
  else {
    solution <- matrix(est.pars[model.set.final$index.matrix], 
                       dim(model.set.final$index.matrix))
    solution.se <- matrix(0, dim(solution)[1], dim(solution)[1])
    eigval <- NULL
    eigvect <- NULL
  }
  if ((any(solution == lb, na.rm = TRUE) || any(solution == 
                                                ub, na.rm = TRUE)) && (lb != 0 || ub != 100)) {
    bound.hit <- TRUE
  }
  rownames(solution) <- rownames(solution.se) <- state.names
  colnames(solution) <- colnames(solution.se) <- state.names
  if (is.character(node.states)) {
    if (node.states == "marginal" || node.states == "scaled") {
      colnames(lik.anc$lik.anc.states) <- state.names
    }
  }
  obj = list(loglik = loglik, AIC = -2 * loglik + 2 * model.set.final$np, 
             AICc = -2 * loglik + (2 * model.set.final$np * (nb.tip/(nb.tip - 
                                                                       model.set.final$np - 1))), ntraits = 1, solution = solution, 
             solution.se = solution.se, index.mat = model.set.final$index.matrix, 
             lewis.asc.bias = lewis.asc.bias, opts = opts, data = data, 
             phy = phy, states = lik.anc$lik.anc.states, tip.states = tip.states, 
             iterations = out$iterations, eigval = eigval, eigvect = eigvect, 
             bound.hit = bound.hit, 
             lower = lb, upper = ub, par.vec = est.pars, root.p = root.p)
  if (!is.null(matching$message.data)) {
    obj$message.data <- matching$message.data
    obj$data <- matching$data
  }
  if (!is.null(matching$message.tree)) {
    obj$message.tree <- matching$message.tree
    obj$data <- matching$data
  }
  class(obj) <- "raydisc"
  return(obj)
}
