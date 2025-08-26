

#' Title Calculate likilehood for Sinba Model (coorelated evolution of two characters)

p 
phy
liks
rate
root.p

#Qindv.zero, Pindv, edge.pattern.focal, cub.method = "pcubature"

likelihood_hmm_multi <- function(p, phy, liks, rate, root.p) {
  
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  anc <- unique(phy$edge[, 1])
  
  # if (is.null(rate)) {
  #   Q <- Qindv.zero
  # } else {
  #   if (any(is.nan(p)) || any(is.infinite(p))) {
  #     return(1e+06)
  #   }
  #   
  #   Q <- mapply(function(x, y) {
  #     x[] <- c(p, 0)[y]
  #     diag(x) <- -rowSums(x)
  #     x
  #   }, x = Qindv.zero, y = rate, SIMPLIFY = FALSE)
  # }
  
  if (any(is.nan(p)) || any(is.infinite(p))) {
    return(1e+06)
  }
  
  # Q <- initQ(c(1, 2), c(1, 1), diag.as = 0)
  # Q <- mapply(function(x, y) {
  #   x[] <- c(p, 0)[y]
  #   diag(x) <- -rowSums(x)
  #   x
  # }, x = Q, y = rate[[1]], SIMPLIFY = FALSE)
  # Q[] <- c(p, 0)[rate[[1]]]
  # diag(Q) <- -rowSums(Q)
  
  #i <- 6
  for (i in seq(from = 1, length.out = nb.node)) {
    focal <- anc[i]
    desRows <- which(phy$edge[, 1] == focal)
    desNodes <- phy$edge[desRows, 2]
    v <- 1
    # desIndex <-1
    for (desIndex in sequence(length(desRows))) {
      # v <- v * expm::expm(Q * phy$edge.length[desRows[desIndex]], method = c("Ward77")) %*% liks[desNodes[desIndex],]
      LL <- NULL
      
      # Tau <- phy$edge.length[desRows[desIndex]]
      # > patternList()
      # [,1]   [,2]   [,3]   [,4]   [,5]     [,6]     [,7] [,8] [,9] [,10] [,11]
      # pat.names "P1Q1" "P2Q2" "Q1Qc" "Q2Qc" "P1Q1Qc" "P2Q2Qc" "Qc" "Qc" "Q1" "Q2"  "P00"
      # pat.vals  "*1"   "*2"   "1*2"  "2*1"  "*1*2"   "*2*1"   "21" "12" "1"  "2"   "0"
      # pat.ids   "1"    "2"    "3"    "4"    "5"      "5"      "6"  "6"  "7"  "8"   "9"
      
      # plot(scen$tree)
      # nodelabels()
      # edgelabels()
      # edgelabels(scen$edge.events[[1]], c(1:nrow(tree$edge)), cex = 0.7)
      
      
      if (edge.pattern.focal[desRows[desIndex]] == 1) { #
        Tau <- phy$edge.length[desRows[desIndex]]
        # LL <- Pindv$P1 %*% ((1/Tau)*integrateOneBirth(Tau=Tau, Q$Q1, cub.method=cub.method)) %*% liks[desNodes[desIndex],]
        LL <- Pindv$P1 %*% ((1 / Tau) * integrateOneBirth_C(Tau = Tau, Q$Q1, cub.method = cub.method)) %*% liks[desNodes[desIndex], ]
      } else if (edge.pattern.focal[desRows[desIndex]] == 2) { #
        Tau <- phy$edge.length[desRows[desIndex]]
        # LL <- Pindv$P2 %*% ((1/Tau)*integrateOneBirth(Tau=Tau, Q$Q2, cub.method=cub.method)) %*% liks[desNodes[desIndex],]
        LL <- Pindv$P2 %*% ((1 / Tau) * integrateOneBirth_C(Tau = Tau, Q$Q2, cub.method = cub.method)) %*% liks[desNodes[desIndex], ]
      } else if (edge.pattern.focal[desRows[desIndex]] == 3) { #
        Tau <- phy$edge.length[desRows[desIndex]]
        # LL <- 1/Tau*integrateScen3(Tau=Tau, Q1=Q$Q1, Q2=Q$Qcom, cub.method=cub.method) %*% liks[desNodes[desIndex],]
        LL <- 1 / Tau * integrateScen3_C(Tau = Tau, Q1 = Q$Q1, Q2 = Q$Qcom, cub.method = cub.method) %*% liks[desNodes[desIndex], ]
      } else if (edge.pattern.focal[desRows[desIndex]] == 4) { #
        Tau <- phy$edge.length[desRows[desIndex]]
        # LL <- 1/Tau*integrateScen3(Tau=Tau, Q1=Q$Q2, Q2=Q$Qcom, cub.method=cub.method) %*% liks[desNodes[desIndex],]
        LL <- 1 / Tau * integrateScen3_C(Tau = Tau, Q1 = Q$Q2, Q2 = Q$Qcom, cub.method = cub.method) %*% liks[desNodes[desIndex], ]
      } else if (edge.pattern.focal[desRows[desIndex]] == 5) { #
        Tau <- phy$edge.length[desRows[desIndex]]
        # MM <- (Pindv$P1 %*% (1/Tau^2*integrate2_Exp(Q1=Q$Q1, Q2=Q$Qcom, Tau=Tau, cub.method=cub.method))) +
        #               (Pindv$P2 %*% (1/Tau^2*integrate2_Exp(Q1=Q$Q2, Q2=Q$Qcom, Tau=Tau, cub.method=cub.method)))
        MM <- (Pindv$P1 %*% (1 / Tau^2 * integrate2_Exp_C(Q1 = Q$Q1, Q2 = Q$Qcom, Tau = Tau, cub.method = cub.method))) +
          (Pindv$P2 %*% (1 / Tau^2 * integrate2_Exp_C(Q1 = Q$Q2, Q2 = Q$Qcom, Tau = Tau, cub.method = cub.method)))
        LL <- MM %*% liks[desNodes[desIndex], ]
      } else if (edge.pattern.focal[desRows[desIndex]] == 6) { #
        # LL <- expm::expm(Q$Qcom * phy$edge.length[desRows[desIndex]], method = c("Ward77")) %*% liks[desNodes[desIndex],]
        LL <- ExpMat_C(phy$edge.length[desRows[desIndex]], Q$Qcom) %*% liks[desNodes[desIndex], ]
      } else if (edge.pattern.focal[desRows[desIndex]] == 7) { #
        # LL <- expm::expm(Q$Q1 * phy$edge.length[desRows[desIndex]], method = c("Ward77")) %*% liks[desNodes[desIndex],]
        LL <- ExpMat_C(phy$edge.length[desRows[desIndex]], Q$Q1) %*% liks[desNodes[desIndex], ]
      } else if (edge.pattern.focal[desRows[desIndex]] == 8) { #
        # LL <- expm::expm(Q$Q2 * phy$edge.length[desRows[desIndex]], method = c("Ward77")) %*% liks[desNodes[desIndex],]
        LL <- ExpMat_C(phy$edge.length[desRows[desIndex]], Q$Q2) %*% liks[desNodes[desIndex], ]
      } else if (edge.pattern.focal[desRows[desIndex]] == 9) { #
        LL <- Pindv$P0 %*% liks[desNodes[desIndex], ]
      }
      
      v <- v * LL
    }
    comp[focal] <- sum(v)
    liks[focal, ] <- v / comp[focal]
    #liks[focal, ] <- v
  }
  
  root <- nb.tip + 1L
  
  if (is.na(sum(log(comp[-TIPS])))) {
    return(1e+06)
  } else {
    if (is.character(root.p)) {
      if (root.p == "maddfitz") {
        root.p <- liks[root, ] / sum(liks[root, ])
      }
    }
    
    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p) + log(liks[root, ])))))
  }
  
  if (is.infinite(loglik)) {
    return(1e+06)
  }
  
  return(loglik)
}