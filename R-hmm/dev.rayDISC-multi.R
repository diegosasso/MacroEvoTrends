library(Rcpp)
library(RcppArmadillo)
# Compile the C++ code
Rcpp::sourceCpp("R-hmm/src/mat.cpp")

source('R-hmm/Lewis.R')

# p <- log(out$solution)
# # liks = model.set.final$liks
# char_env[['char1']]
# liks <- char_env
# Q = model.set.final$Q
# rate = model.set.final$rate
# 
# root.p = root.p
# lewis.asc.bias = lewis.asc.bias

dev.raydisc_multi <- function (p, phy, char_env, char_invar_env, Q, rate, root.p, lewis.asc.bias) 
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
  # Check row sums ≈ 0
  row_sums <- rowSums(Q)
  tol = 1e-10
  if (any(abs(row_sums) > tol)) {
    print(Q)
    stop("Invalid Q: row sums are not approximately zero.")
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
  # #------- Traverse
  # comp <- numeric(nb.tip + nb.node)
  # i=1
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
  #
  #--------- ROOT
  root <- nb.tip + 1L
  # Determine root.p based on user input or defaults
  nstates <- ncol(Q)
  if (is.character(root.p)) {
    root.p <- tolower(root.p)
    
    if (root.p == "flat") {
      root.p <- rep(1/nstates, nstates)
      
    } else if (root.p == "yang") {
      # Use stationary distribution from Q
      root.p <- MASS::Null(Q)
      root.p <- c(root.p / sum(root.p))
      
    } else if (root.p == "madfitz") {
      # Use scaled likelihoods at the root
      root.p <- liks[root, ] / sum(liks[root, ])
      
    } else {
      stop("Unknown root.p option. Use 'flat', 'yang', 'madfitz', or a numeric vector.")
    }
  } 
  else if (is.numeric(root.p)) {
    # Normalize custom root.p to sum to 1
    root.p <- root.p / sum(root.p)
  } 
  else {
    stop("Invalid root.p: must be NULL, 'flat', 'yang', 'madfitz', or numeric.")
  }
  #---------- END ROOT
  #
  #---------- Traverse
  # Prepare results storage
  char_names <- ls(char_env)
  Nchar <- length(char_names)
  
  # Compensation factors per character
  comp_list <- vector("list", Nchar)
  names(comp_list) <- char_names
  loglik_total <- 0
  # Likelihoods updated per character (environment already has tip likelihoods)
  # ci=1
  for (ci in seq_along(char_names)) {
    char <- char_names[ci]
    liks <- char_env[[char]]  # matrix [tips × states]
    
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
          # v <- v * expm(Q * phy$edge.length[desRows[desIndex]], method = c("Ward77")) %*% liks[desNodes[desIndex],]
          v <- v * ExpMat_C(phy$edge.length[desRows[desIndex]], Q) %*% liks[desNodes[desIndex],]
        }
        
        # Store compensation and normalized likelihood
        comp[focal]   <- sum(v)
        liks[focal, ] <- v / comp[focal]
    }
    
    # Compute per-character log-likelihood
    loglik_char <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root, ]))))
    
    # Handle numerical underflow or bad values
    if (is.infinite(loglik_char) || is.na(loglik_char)) {
      return(1e6)  # Large penalty if invalid likelihood
    }
    
    # Accumulate total log-likelihood
    loglik_total <- loglik_total + loglik_char
    
    # Save updated likelihoods back into environment
    # char_env[[char]] <- liks
    # comp_list[[ci]]  <- comp
  }
 #  #--------- ROOT
 #  #
 #  root <- nb.tip + 1L
 #  # if (is.na(sum(log(comp[-TIPS])))) {
 #  comp.isna <- lapply(comp_list, function(x) is.na(sum(log(x[-TIPS]))))
 #  comp.isna <- unlist(comp.isna)
 #  if (any(comp.isna)) {
 #    return(1e+06)
 #  } 
 #  
 # # Determine root.p based on user input or defaults
 # nstates <- ncol(Q)
 # if (is.character(root.p)) {
 #    root.p <- tolower(root.p)
 #    
 #    if (root.p == "flat") {
 #      root.p <- rep(1/nstates, nstates)
 #      
 #    } else if (root.p == "yang") {
 #      # Use stationary distribution from Q
 #      root.p <- MASS::Null(Q)
 #      root.p <- c(root.p / sum(root.p))
 #      
 #    } else if (root.p == "madfitz") {
 #      # Use scaled likelihoods at the root
 #      root.p <- liks[root, ] / sum(liks[root, ])
 #      
 #    } else {
 #      stop("Unknown root.p option. Use 'flat', 'yang', 'madfitz', or a numeric vector.")
 #    }
 # } 
 # else if (is.numeric(root.p)) {
 #    # Normalize custom root.p to sum to 1
 #    root.p <- root.p / sum(root.p)
 # } 
 # else {
 #    stop("Invalid root.p: must be NULL, 'flat', 'yang', 'madfitz', or numeric.")
 # }
 
 
 #---------------------- ORIGINAL ROOT
  # if (is.na(sum(log(comp[-TIPS])))) {
  #   return(1e+06)
  # }
  # else {
  #   equil.root <- NULL
  #   for (i in 1:ncol(Q)) {
  #     posrows <- which(Q[, i] >= 0)
  #     rowsum <- sum(Q[posrows, i])
  #     poscols <- which(Q[i, ] >= 0)
  #     colsum <- sum(Q[i, poscols])
  #     equil.root <- c(equil.root, rowsum/(rowsum + colsum))
  #   }
  #   if (is.null(root.p)) {
  #     flat.root = equil.root
  #     k.rates <- 1/length(which(!is.na(equil.root)))
  #     flat.root[!is.na(flat.root)] = k.rates
  #     flat.root[is.na(flat.root)] = 0
  #     root.p <- flat.root
  #     loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root, ]))))
  #   }
  #   else {
  #     if (is.character(root.p)) {
  #       # stationary root
  #       if (root.p == "yang") {
  #         root.p <- MASS::Null(Q)
  #         #expm(Q*100)
  #         root.p <- c(root.p/sum(root.p))
  #         loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root, ]))))
  #         if (is.infinite(loglik)) {
  #           return(1e+06)
  #         }
  #       }
  #       # maddfitz
  #       else {
  #         root.p = liks[root, ]/sum(liks[root, ])
  #         loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root, ]))))
  #       }
  #     }
  #     # custom numeric root
  #     else {
  #       loglik <- sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root, ]))))
  #       if (is.infinite(loglik)) {
  #         return(1e+06)
  #       }
  #     }
  #   }
  # }

  
  #---- lewis.asc.bias
  if (lewis.asc.bias == TRUE) {
    # dummy.liks.vec <- numeric(dim(Q)[1])
    # for (state.index in 1:dim(Q)[1]) {
    #   dummy.liks.vec[state.index] <- corHMM:::CalculateLewisLikelihood(p = p.new, 
    #                                                           phy = phy, liks = liks, Q = Q, rate = rate, 
    #                                                           root.p = root.p, state.num = state.index)
    # }
    # loglik <- loglik - log(sum(root.p * (1 - exp(dummy.liks.vec))))
    #----------------
    char_names_invar <- ls(char_invar_env)
    #print(char_names_invar)
    dummy.liks.vec <- numeric(length(char_names_invar))
    # ci=1
    for (ci in seq_along(char_names_invar)) {
      char <- char_names_invar[ci]
      #print(char)
      liks.invar <- char_invar_env[[char]]
      dummy.liks.vec[ci] <- ProbInvarLewis(phy = phy, liks = liks.invar, Q = Q, root.p = root.p)
    }
    # dummy.liks.vec <- numeric(dim(Q)[1])
    # for (state.index in 1:dim(Q)[1]) {
    #   dummy.liks.vec[state.index] <- ProbInvarLewis(phy = phy, liks = liks, Q = Q, root.p = root.p, state.num = state.index)
    # }
    LewisCorr <- log(1-sum(exp(dummy.liks.vec))) * Nchar
    loglik_total <- loglik_total - LewisCorr
  }
  #
  -loglik_total
}

# log(0.003)
# dd <- c(-1.501552, -9.453301)
# log(sum(1-exp(dd)))
# log( 50 * (1-sum(exp(dd))) )
# log(  (1-sum(exp(dd))) ) +log(50)
