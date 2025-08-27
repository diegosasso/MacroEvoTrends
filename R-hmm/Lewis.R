ProbInvarLewis <- function (phy, liks, Q, root.p) 
{
  #p.new <- p
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  anc <- unique(phy$edge[, 1])
  # if (is.null(rate)) {
  #   Q = Q
  # }
  # else {
  #   if (any(is.nan(p.new)) || any(is.infinite(p.new))) 
  #     return(1e+06)
  #   Q[] <- c(p.new, 0)[rate]
  #   diag(Q) <- -rowSums(Q)
  # }
  liks.dummy <- liks
  #liks.dummy[TIPS, ] = 0
  #liks.dummy[ , ] <- 0
  #liks.dummy[TIPS, state.num] = 1
  comp.dummy <- comp
  for (i in seq(from = 1, length.out = nb.node)) {
    focal <- anc[i]
    desRows <- which(phy$edge[, 1] == focal)
    desNodes <- phy$edge[desRows, 2]
    v.dummy <- 1
    for (desIndex in sequence(length(desRows))) {
      # v.dummy <- v.dummy * expm(Q * phy$edge.length[desRows[desIndex]], method = c("Ward77")) %*% liks.dummy[desNodes[desIndex], ]
      v.dummy <- v.dummy * ExpMat_C(phy$edge.length[desRows[desIndex]], Q) %*% liks.dummy[desNodes[desIndex], ]
    }
    comp.dummy[focal] <- sum(v.dummy)
    liks.dummy[focal, ] <- v.dummy/comp.dummy[focal]
  }
  root <- nb.tip + 1L
  if (is.na(sum(log(comp[-TIPS])))) {
    return(1e+06)
  }
  loglik <- (sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(root.p) + log(liks.dummy[root, ])))))
  if (is.infinite(loglik)) {
    return(1e+06)
  }
  loglik
}
