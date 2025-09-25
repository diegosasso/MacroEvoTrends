# This is a modified version of corHMM funciton that returns maddfitz root values
#
ancRECON_ROOT <- function (phy, data, p, method = c("joint", "marginal", "scaled"), 
          rate.cat, ntraits = NULL, rate.mat = NULL, model = "ARD", 
          root.p = NULL, get.likelihood = FALSE, get.tip.states = FALSE, 
          collapse = TRUE) 
{
  if (!is.null(root.p)) {
    if (!is.character(root.p)) {
      root.p <- root.p/sum(root.p)
    }
  }
  root.p_input <- root.p
  if (is.null(rate.cat)) {
    rate.cat <- 1
  }
  input.data <- data
  corData <- corHMM:::corProcessData(data, collapse = collapse)
  data <- corData$corData
  matching <- corHMM:::match.tree.data(phy, data)
  data <- matching$data
  phy <- matching$phy
  phy$edge.length[phy$edge.length <= 1e-05] = 1e-05
  data.sort <- data.frame(data[, 2], data[, 2], row.names = data[, 
                                                                 1])
  data.sort <- data.sort[phy$tip.label, ]
  levels <- levels(as.factor(data.sort[, 1]))
  obj <- NULL
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  k <- ntraits <- length(corData$ObservedTraits)
  drop.states = NULL
  if (is.null(rate.mat)) {
    model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy = phy, 
                                               data = input.data, rate.cat = rate.cat, ntraits = ntraits, 
                                               model = model, rate.mat = rate.mat, collapse = collapse)
    rate.mat <- model.set.final$index.matrix
    rate <- model.set.final$rate
  }
  else {
    model.set.final <- rate.cat.set.corHMM.JDB(phy = phy, 
                                               data = input.data, rate.cat = rate.cat, ntraits = ntraits, 
                                               model = model, rate.mat = rate.mat, collapse = collapse)
    rate <- rate.mat
    col.sums <- which(colSums(rate.mat, na.rm = TRUE) == 
                        0)
    row.sums <- which(rowSums(rate.mat, na.rm = TRUE) == 
                        0)
    drop.states <- col.sums[which(col.sums == row.sums)]
    rate[is.na(rate)] <- max(rate, na.rm = TRUE) + 1
  }
  x <- data.sort[, 1]
  TIPS <- 1:nb.tip
  tranQ <- Q <- model.set.final$Q
  liks <- model.set.final$liks
  if (length(drop.states > 0)) {
    liks[, drop.states] <- 0
  }
  p[p == 0] = exp(-21)
  Q[] <- c(p, 0)[rate]
  diag(Q) <- -rowSums(Q)
  phy <- reorder(phy, "pruningwise")
  TIPS <- 1:nb.tip
  anc <- unique(phy$edge[, 1])
  if (method == "joint") {
    if (!is.null(phy$node.label)) {
      tip.state.vector <- rep(NA, Ntip(phy))
      known.state.vector <- phy$node.label
      known.state.vector <- c(tip.state.vector, known.state.vector)
    }
    else {
      tip.state.vector <- rep(NA, Ntip(phy))
      known.state.vector <- rep(NA, Nnode(phy))
      known.state.vector <- c(tip.state.vector, known.state.vector)
    }
    lik.states <- numeric(nb.tip + nb.node)
    pupko.L <- matrix(NA, nrow = nb.tip + nb.node, ncol(liks))
    pupko.C <- matrix(NA, nrow = nb.tip + nb.node, ncol(liks))
    for (i in seq(from = 1, length.out = nb.node)) {
      focal <- anc[i]
      desRows <- which(phy$edge[, 1] == focal)
      desNodes <- phy$edge[desRows, 2]
      for (desIndex in sequence(length(desRows))) {
        if (any(desNodes[desIndex] == phy$edge[, 1]) == 
            FALSE) {
          v <- c(rep(1, k * rate.cat))
          Pij <- expm::expm(Q * phy$edge.length[desRows[desIndex]], 
                      method = c("Ward77"))
          v <- v * liks[desNodes[desIndex], ]
          L <- Pij %*% v
          if (is.na(known.state.vector[focal])) {
            pupko.L[desNodes[desIndex], ] <- L
            pupko.C[desNodes[desIndex], ] <- which.is.max(L == 
                                                            max(L))
          }
          else {
            pupko.L[desNodes[desIndex], ] <- L[known.state.vector[focal], 
            ]
            pupko.C[desNodes[desIndex], ] <- known.state.vector[focal]
          }
        }
      }
      tz <- phy$edge.length[which(phy$edge[, 2] == focal)]
      if (length(tz) == 0) {
        root.state = 1
        for (desIndex in sequence(length(desRows))) {
          root.state <- root.state * pupko.L[desNodes[desIndex], 
          ]
        }
        if (is.na(known.state.vector[focal])) {
          equil.root <- NULL
          for (i in 1:ncol(Q)) {
            posrows <- which(Q[, i] >= 0)
            rowsum <- sum(Q[posrows, i])
            poscols <- which(Q[i, ] >= 0)
            colsum <- sum(Q[i, poscols])
            equil.root <- c(equil.root, rowsum/(rowsum + 
                                                  colsum))
          }
          if (is.null(root.p)) {
            if (is.na(known.state.vector[focal])) {
              flat.root = equil.root
              k.rates <- 1/length(which(!is.na(equil.root)))
              flat.root[!is.na(flat.root)] = k.rates
              flat.root[is.na(flat.root)] = 0
              root.p <- flat.root
              pupko.L[focal, ] <- root.state
            }
            else {
              root.p <- rep(0, dim(Q)[2])
              root.p[known.state.vector[focal]] <- 1
              pupko.L[focal, ] <- root.state
            }
          }
          else {
            if (is.character(root.p)) {
              if (root.p == "yang") {
                root.p <- Null(Q)
                root.p <- c(root.p/sum(root.p))
                pupko.L[focal, ] <- root.state
              }
              else {
                root.p <- (root.state/sum(root.state))[]
                pupko.L[focal, ] <- root.state
              }
            }
            else {
              root.p <- root.p
              pupko.L[focal, ] <- root.state
            }
          }
        }
        else {
          root.p = rep(0, dim(Q)[1])
          root.p[known.state.vector[focal]] <- 1
          pupko.L[focal, ] <- root.state
        }
      }
      else {
        Pij <- expm::expm(Q * tz, method = c("Ward77"))
        v <- c(rep(1, k * rate.cat))
        if (is.na(known.state.vector[focal])) {
          for (desIndex in sequence(length(desRows))) {
            v <- v * pupko.L[desNodes[desIndex], ]
          }
          focalRow <- which(phy$edge[, 2] == focal)
          motherRow <- which(phy$edge[, 1] == phy$edge[focalRow, 
                                                       1])
          motherNode <- phy$edge[focalRow, 1]
          if (is.na(known.state.vector[motherNode])) {
            for (row.index in 1:dim(Pij)[1]) {
              L <- Pij[row.index, ] * v
              pupko.L[focal, row.index] <- max(L)
              pupko.C[focal, row.index] <- which.is.max(L)
            }
          }
          else {
            L <- Pij[known.state.vector[motherNode], 
            ] * v
            pupko.L[focal, ] <- L
            pupko.C[focal, ] <- which.is.max(L)
          }
        }
        else {
          for (desIndex in sequence(length(desRows))) {
            v <- v * pupko.L[desNodes[desIndex], ]
          }
          focalRow <- which(phy$edge[, 2] == focal)
          motherRow <- which(phy$edge[, 1] == phy$edge[focalRow, 
                                                       1])
          motherNode <- phy$edge[focalRow, 1]
          if (is.na(known.state.vector[motherNode])) {
            for (row.index in 1:dim(Pij)[1]) {
              L <- Pij[row.index, ] * v
              pupko.L[focal, row.index] <- L[known.state.vector[focal]]
              pupko.C[focal, row.index] <- known.state.vector[focal]
            }
          }
          else {
            L <- Pij[known.state.vector[motherNode], 
            ] * v
            pupko.L[focal, ] <- L[known.state.vector[focal]]
            pupko.C[focal, ] <- known.state.vector[focal]
          }
        }
        if (sum(pupko.L[focal, ]) < 1e-200) {
          cat("Kicking in arbitrary precision package Rmpfr due to very low probabilities.\n")
          pupko.L <- mpfr(pupko.L, 15)
        }
      }
    }
    root <- nb.tip + 1L
    if (get.likelihood == TRUE) {
      loglik <- log(sum(exp(log(root.p) + log(pupko.L[root, 
      ]))))
      return(as.numeric(loglik))
    }
    else {
      root <- nb.tip + 1L
      if (is.na(known.state.vector[root])) {
        pupko.L[root, ] <- log(root.p) + log(pupko.L[root, 
        ])
        lik.states[root] <- which(pupko.L[root, ] == 
                                    max(pupko.L[root, ]))[1]
      }
      else {
        lik.states[root] <- known.state.vector[root]
      }
      N <- dim(phy$edge)[1]
      for (i in N:1) {
        anc <- phy$edge[i, 1]
        des <- phy$edge[i, 2]
        lik.states[des] <- pupko.C[des, lik.states[anc]]
      }
      obj$lik.tip.states <- lik.states[TIPS]
      obj$lik.anc.states <- lik.states[-TIPS]
      obj$info.anc.states <- NULL
      return(obj)
    }
  }
  if (method == "marginal") {
    liks.down <- liks
    tranQ <- t(Q)
    comp <- matrix(0, nb.tip + nb.node, ncol(liks))
    for (i in seq(from = 1, length.out = nb.node)) {
      focal <- anc[i]
      desRows <- which(phy$edge[, 1] == focal)
      desNodes <- phy$edge[desRows, 2]
      v <- 1
      for (desIndex in sequence(length(desRows))) {
        v <- v * expm::expm(Q * phy$edge.length[desRows[desIndex]], 
                      method = c("Ward77")) %*% liks.down[desNodes[desIndex], 
                      ]
      }
      if (!is.null(phy$node.label)) {
        if (!is.na(phy$node.label[focal - nb.tip])) {
          fixer.tmp = numeric(dim(Q)[2]/rate.cat)
          fixer.tmp[phy$node.label[focal - nb.tip]] = 1
          fixer = rep(fixer.tmp, rate.cat)
          v <- v * fixer
        }
      }
      comp[focal] <- sum(v)
      liks.down[focal, ] <- v/comp[focal]
    }
    root <- nb.tip + 1L
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
      liks.down[root, ] <- flat.root * liks.down[root, 
      ]
      liks.down[root, ] <- liks.down[root, ]/sum(liks.down[root, 
      ])
      root.p = flat.root
    }
    else {
      if (is.character(root.p)) {
        if (root.p == "yang") {
          root.p <- Null(Q)
          root.p <- c(root.p/sum(root.p))
          liks.down[root, ] <- liks.down[root, ] * root.p
          liks.down[root, ] <- liks.down[root, ]/sum(liks.down[root, 
          ])
        }
        else {
          root.p = liks.down[root, ]/sum(liks.down[root, 
          ])
          liks.down[root, ] <- root.p * liks.down[root, 
          ]
          liks.down[root, ] <- liks.down[root, ]/sum(liks.down[root, 
          ])
        }
      }
      else {
        liks.down[root, ] <- root.p * liks.down[root, 
        ]
        liks.down[root, ] <- liks.down[root, ]/sum(liks.down[root, 
        ])
      }
    }
    obsRoot <- root.p
    root.p <- vector("numeric", dim(liks)[2])
    root.p[] <- obsRoot
    liks.up <- liks
    states <- apply(liks, 1, which.max)
    N <- dim(phy$edge)[1]
    comp <- numeric(nb.tip + nb.node)
    for (i in length(anc):1) {
      focal <- anc[i]
      if (!focal == root) {
        focalRow <- which(phy$edge[, 2] == focal)
        motherRow <- which(phy$edge[, 1] == phy$edge[focalRow, 
                                                     1])
        motherNode <- phy$edge[focalRow, 1]
        desNodes <- phy$edge[motherRow, 2]
        sisterNodes <- desNodes[(which(!desNodes == 
                                         focal))]
        sisterRows <- which(phy$edge[, 2] %in% sisterNodes == 
                              TRUE)
        if (motherNode != root) {
          v <- expm(tranQ * phy$edge.length[which(phy$edge[, 
                                                           2] == motherNode)], method = c("Ward77")) %*% 
            liks.up[motherNode, ]
          if (!is.null(phy$node.label)) {
            if (!is.na(phy$node.label[motherNode - nb.tip])) {
              fixer.tmp <- numeric(dim(Q)[2]/rate.cat)
              fixer.tmp[phy$node.label[motherNode - 
                                         nb.tip]] <- 1
              fixer <- rep(fixer.tmp, rate.cat)
              v <- v * fixer
            }
          }
        }
        else {
          v <- root.p
        }
        for (sisterIndex in sequence(length(sisterRows))) {
          v <- v * expm(Q * phy$edge.length[sisterRows[sisterIndex]], 
                        method = c("Ward77")) %*% liks.down[sisterNodes[sisterIndex], 
                        ]
        }
        comp[focal] <- sum(v)
        liks.up[focal, ] <- v/comp[focal]
      }
    }
    liks.final <- liks
    comp <- numeric(nb.tip + nb.node)
    for (i in seq(from = 1, length.out = nb.node - 1)) {
      focal <- anc[i]
      focalRows <- which(phy$edge[, 2] == focal)
      v <- liks.down[focal, ] * (expm(tranQ * phy$edge.length[focalRows], 
                                      method = c("Ward77")) %*% liks.up[focal, ])
      comp[focal] <- sum(v)
      liks.final[focal, ] <- v/comp[focal]
    }
    if (get.tip.states == TRUE) {
      liks.final[TIPS, ] <- GetTipStateBruteForce(p = p, 
                                                  phy = phy, data = input.data, rate.mat = rate.mat, 
                                                  rate.cat = rate.cat, ntraits = ntraits, model = model, 
                                                  root.p = root.p_input, collapse = collapse)
    }
    else {
      liks.final[TIPS, ] <- liks.down[TIPS, ]
    }
    liks.final[root, ] <- liks.down[root, ]
    if (get.likelihood == TRUE) {
      loglik <- as.numeric(log(liks[root, lik.states[root]]))
      return(loglik)
    }
    else {
      obj$lik.tip.states <- liks.final[TIPS, ]
      obj$lik.anc.states <- liks.final[-TIPS, ]
      obj$info.anc.states <- getInfoPerNode(obj$lik.anc.states, 
                                            Q)
      return(obj)
    }
  }
  if (method == "scaled") {
    comp <- matrix(0, nb.tip + nb.node, ncol(liks))
    root <- nb.tip + 1L
    for (i in seq(from = 1, length.out = nb.node)) {
      focal <- anc[i]
      desRows <- which(phy$edge[, 1] == focal)
      desNodes <- phy$edge[desRows, 2]
      v <- 1
      for (desIndex in sequence(length(desRows))) {
        v <- v * expm::expm(Q * phy$edge.length[desRows[desIndex]], 
                      method = c("Ward77")) %*% liks[desNodes[desIndex], 
                      ]
      }
      comp[focal] <- sum(v)
      liks[focal, ] <- v/comp[focal]
    }
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
      liks[root, ] <- flat.root * liks[root, ]
      liks[root, ] <- liks[root, ]/sum(liks[root, ])
    }
    else {
      if (is.character(root.p)) {
        if (root.p == "yang") {
          root.p <- Null(Q)
          root.p <- c(root.p/sum(root.p))
          liks[root, ] <- root.p * liks[root, ]
          liks[root, ] <- liks[root, ]/sum(liks[root, 
          ])
        }
        else {
          root.p = liks[root, ]/sum(liks[root, ])
          liks[root, ] <- root.p * liks[root, ]
          liks[root, ] <- liks[root, ]/sum(liks[root, 
          ])
        }
      }
      else {
        liks[root, ] <- root.p * liks[root, ]
        liks[root, ] <- liks[root, ]/sum(liks[root, 
        ])
      }
    }
    obj$lik.tip.states <- liks[TIPS, ]
    obj$lik.anc.states <- liks[-TIPS, ]
    obj$info.anc.states <- NULL
    obj$root <- root.p
    return(obj)
  }
}
