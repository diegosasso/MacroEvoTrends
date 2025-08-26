p <- log(out$solution)
# liks = model.set.final$liks
char_env[['char1']]
liks <- char_env
Q = model.set.final$Q
rate = model.set.final$rate

root.p = root.p 
lewis.asc.bias = lewis.asc.bias

dev.raydisc_multi <- function (p, phy, char_env, Q, rate, root.p, lewis.asc.bias) 
{
  p.new <- exp(p)
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  # comp <- numeric(nb.tip + nb.node)
  anc <- unique(phy$edge[, 1])
  if (is.null(rate)) {
    Q = Q
  }else {
    if (any(is.nan(p.new)) || any(is.infinite(p.new))) 
      return(1e+06)
    Q[] <- c(p.new, 0)[rate]
    diag(Q) <- -rowSums(Q)
  }
  # #------- Traverse
  # for (i in seq(from = 1, length.out = nb.node)) {
  #   focal <- anc[i]
  #   desRows <- which(phy$edge[, 1] == focal)
  #   desNodes <- phy$edge[desRows, 2]
  #   v <- 1
  #   for (desIndex in sequence(length(desRows))) {
  #     v <- v * expm(Q * phy$edge.length[desRows[desIndex]], 
  #                   method = c("Ward77")) %*% liks[desNodes[desIndex], 
  #                   ]
  #   }
  #   comp[focal] <- sum(v)
  #   liks[focal, ] <- v/comp[focal]
  # }
  # #-------
  #------- Traverse
  comp <- numeric(nb.tip + nb.node)
  i=1
  for (i in seq(from = 1, length.out = nb.node)) {
    focal <- anc[i]
    desRows <- which(phy$edge[, 1] == focal)
    desNodes <- phy$edge[desRows, 2]
    v <- 1
    for (desIndex in sequence(length(desRows))) {
      v <- v * expm(Q * phy$edge.length[desRows[desIndex]], 
                    method = c("Ward77")) %*% liks[desNodes[desIndex], 
                    ]
    }
    comp[focal] <- sum(v)
    liks[focal, ] <- v/comp[focal]
  }
  #----------
  #----------
  # Prepare results storage
  char_names <- ls(char_env)
  Nchar <- length(char_names)
  
  # Compensation factors per character
  comp_list <- vector("list", Nchar)
  names(comp_list) <- char_names
  
  # Likelihoods updated per character (environment already has tip likelihoods)
  # ci=1
  for (ci in seq_along(char_names)) {
    char <- char_names[ci]
    liks <- char_env[[char]]  # matrix [tips × states]
    #Q    <- Q_list[[ci]]      # substitution matrix for this character
    
    # Preallocate comp for this character
    comp <- numeric(nb.tip + nb.node)
    
    # Traverse internal nodes
    #i=1
    for (i in seq(from = 1, length.out = nb.node)) {
      focal    <- anc[i]
      desRows  <- which(phy$edge[, 1] == focal)
      desNodes <- phy$edge[desRows, 2]
      
      v <- 1
      for (desIndex in sequence(length(desRows))) {
        v <- v * expm(Q * phy$edge.length[desRows[desIndex]], method = c("Ward77")) %*% liks[desNodes[desIndex],]
      }
      
      # Store compensation and normalized likelihood
      comp[focal]   <- sum(v)
      liks[focal, ] <- v / comp[focal]
    }
    
    # Save updated likelihoods back into environment
    char_env[[char]] <- liks
    comp_list[[ci]]  <- comp
  }
  #---------
  #
  root <- nb.tip + 1L
  if (is.na(sum(log(comp[-TIPS])))) {
    return(1e+06)
  }
  else {
    equil.root <- NULL
    for (i in 1:ncol(Q)) {
      posrows <- which(Q[, i] >= 0)
      rowsum <- sum(Q[posrows, i])
      poscols <- which(Q[i, ] >= 0)
      colsum <- sum(Q[i, poscols])
      equil.root <- c(equil.root, rowsum/(rowsum + colsum))
    }
    if (is.null(root.p)) {
      flat.root = equil.root
      k.rates <- 1/length(which(!is.na(equil.root)))
      flat.root[!is.na(flat.root)] = k.rates
      flat.root[is.na(flat.root)] = 0
      root.p <- flat.root
      loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + 
                                                      log(liks[root, ]))))
    }
    else {
      if (is.character(root.p)) {
        # stationary root
        if (root.p == "yang") {
          root.p <- MASS::Null(Q)
          #expm(Q*100)
          root.p <- c(root.p/sum(root.p))
          loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + 
                                                          log(liks[root, ]))))
          if (is.infinite(loglik)) {
            return(1e+06)
          }
        }
        # maddfitz
        else {
          root.p = liks[root, ]/sum(liks[root, ])
          loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + 
                                                          log(liks[root, ]))))
        }
      }
      else {
        loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + 
                                                        log(liks[root, ]))))
        if (is.infinite(loglik)) {
          return(1e+06)
        }
      }
    }
  }
  if (lewis.asc.bias == TRUE) {
    dummy.liks.vec <- numeric(dim(Q)[1])
    for (state.index in 1:dim(Q)[1]) {
      dummy.liks.vec[state.index] <- corHMM:::CalculateLewisLikelihood(p = p.new, 
                                                              phy = phy, liks = liks, Q = Q, rate = rate, 
                                                              root.p = root.p, state.num = state.index)
    }
    loglik <- loglik - log(sum(root.p * (1 - exp(dummy.liks.vec))))
  }
  -loglik
}
