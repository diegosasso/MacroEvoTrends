library(stringr)

is.sequential <- function(x, eps = 1e-8) {
  if (length(x) && isTRUE(abs(x[1] - floor(x[1])) < eps)) {
    all(abs(diff(x) - 1) < eps)
  } else {
    FALSE
  }
}

mk_RatePars <- function(Qindv) {
  uniq.rt <- lapply(Qindv, function(x) c(x)) %>%
    unlist() %>%
    unique()
  uniq.rt <- uniq.rt[uniq.rt > 0]
  
  zero.rate <- max(uniq.rt + 1)
  
  if (!is.sequential(sort(uniq.rt))) stop("Rate paremeters are not sequential")
  
  out <- lapply(Qindv, function(x) {
    x[x == 0] <- zero.rate
    x
  })
  
  return(out)
}


#taxa
make_lik_bins <- function(taxa, Qcom){
  M.nstates <- ncol(Qcom)
  M.ntax <- nrow(taxa)
  
  Ln.bins <- matrix(0, ncol = M.nstates, nrow = M.ntax)
  colnames(Ln.bins) <- colnames(Qcom)
  
  dt <- taxa[,2]
  
  i <- 1
  for (i in 1:length(dt)){
    token <- as.character(dt[i])
    # token <- '1&2&3'
    token <- str_split(token, '&')[[1]] %>% as.numeric()
    Ln.bins[i,token] <- 1
    #Ln.bins[,i][str_detect(dt, token)] <- 1
  }
  
  rownames(Ln.bins) <- rownames(taxa)
  return(Ln.bins)
}

make_lik_bins_all <- function(taxa, Qcom, phy){
  
  Tips.ln <- make_lik_bins(taxa, Qcom)
  rownames(Tips.ln) <- NULL
  out <- rbind(Tips.ln, matrix(0, ncol = ncol(Qcom), nrow = Nnode(phy)) )
  return(out)
}



#' Title Infer corrleated evolution of two characters using Sinba model

tree <- readRDS("tree_test.RDS")
char <- readRDS('RDS/char.rds')

phy <- tree
data <- cbind(rownames(char), char[,1]+1)

Q <- initQ(c(1, 2), c(1, 1), diag.as = 0)
root.p = NULL
#p = NULL
# p=initQ(c(1, 2), c(.1, .1))
p= c(0.1, 0.1)
ip = NULL
lb = 0
ub = 100
diagn=FALSE
optim.method='subplex'
subplex.rounds=1
subplex.n.cores=1L

# edge.pattern=NULL
# Edge.Scenarios
# Pindv=NULL
# cub.method = "pcubature"
# gensa.its=100

fit_hmm_multi <- function(phy, data, root.p = NULL, Q, p = NULL, ip = NULL, lb = 0, ub = 100,
                     diagn=FALSE, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100){
  
  
  
  
  cat("\n/---------------------------------------------------\n")
  cat('Running SinBA model:\n')
  # aa <- c(missing(phy), missing(data))
  # 
  # if (any(aa) &&  missing(Edge.Scenarios))
  #   return(cat('Specify: tree, character, and possibly edge pattern OR edge_scenario object from getEdgeScenarios().\n'))
  # if (any(!aa) &&  !missing(Edge.Scenarios))
  #   return(cat('Specify either (1) tree, character, and possibly edge pattern, OR (2) edge_scenario object from getEdgeScenarios().\n'))
  # 
  # if (!missing(Edge.Scenarios)){
  #   phy <- Edge.Scenarios$tree
  #   data<- Edge.Scenarios$data
  #   
  #   
  #   if(is.null(edge.pattern)){
  #     cat('Edge pattern is not specified. Using the traditional likelihood calculation (no trait births along the edges).\n')
  #     edge.pattern <- rep(6, length(phy$edge.length))
  #   } else if (length(edge.pattern)==1){
  #     cat(paste0('Edge pattern is specified, using edge scenario #', edge.pattern, '.\n'))
  #     edge.pattern <- Edge.Scenarios$edge.patterns[[edge.pattern]]
  #   } else if (length(edge.pattern)>1){
  #     cat('Edge pattern is specified.\n')
  #   }
  #   
  # } else if (missing(Edge.Scenarios)){
  #   
  #   if(is.null(edge.pattern)){
  #     cat('Edge pattern is not specified. Using the traditional likelihood calculation (no trait births along the edges).\n')
  #     edge.pattern <- rep(6, length(phy$edge.length))
  #   } else if (length(edge.pattern)==1){
  #     return(cat('Error: edge pattern is incorrectly specified, please provide it.\n'))
  #   } else if (length(edge.pattern)>1){
  #     cat('Edge pattern is specified.\n')
  #   }
  # }
  
  # # Recode edge patern if neccesary
  # if (is.character(edge.pattern)){
  #   edge.pattern <-edgePatRecode(edge.pattern, from='codes', to='ids')
  # }
  # 
  # #mm <- match(sort(unique(edge.pattern)), as.numeric(patternList()[3,]) )
  # mm <- edgePatRecode(edge.pattern, from='ids', to='codes')
  # cat(paste0('Uniques edge scenarios (codes): ', paste0(unique(mm), collapse = ', '), '.\n' ))
  
  
  # reoder
  #nb.tip <- length(phy$tip.label)
  #TIPS <- 1:nb.tip
  #reorder.phylo
  #ape:::.reorder_ape
  neworder <- reorder(phy, order = "pruningwise", index.only = TRUE)
  
  phy$edge <- phy$edge[neworder, ]
  phy$edge.length <- phy$edge.length[neworder]
  attr(phy, "order") <- "pruningwise"
  #edge.pattern <- edge.pattern[neworder]
  
  if (!all(phy$tip.label %in% data[,1]))
    return('Error: not all taxa in character are present in the tree.\n')
  
  orr <- match(data[,1],phy$tip.label)
  data <- data[orr,]
  # phy <- reorder(phy, order = "pr")
  #reorder(phy, order = "pr", index.only = T)
  
  taxa <- data[, c(1, 2)]
  
  if(is.null(colnames(Q))){
    cat('No row/column names are given in the rate matrix. Naming rows/columns sequentially.\n')
    new.nm <- c(1:ncol(Qindv$Qcom))
    lapply(Qindv, function(x) {colnames(x) <- rownames(x) <-new.nm; x} )
    lapply(Pindv, function(x) {colnames(x) <- rownames(x) <-new.nm; x} )
  }
  
  if (!is.sequential(as.numeric((colnames(Q)))))
    return(cat('Error: Rows/columns in rate matrix are not named sequentially. Rename them.\n'))
  
  #-- checking if states match Q names
  #tk <- c(1,2,3, '5&6', '1&2')
  tk <-taxa[,2]
  tk <- stringr::str_split(tk, '&' ) %>% unlist %>%  unique()
  
  if (!all(tk %in% colnames(Q)))
    return(cat('Error: there are more states in the character than specified by the rate matrix.\n'))
  
  if (any(tk=='0')){return(cat('Error: character states should start with 1 but NOT 0. You need to recode them.\n'))}
  #--
  
  #rate <- mk_RatePars(Qindv)
  #Qindv.zero <- mk_QindvZero(Qindv)
  # !!!
  rate <- mk_RatePars(list(Q))
  
  if (is.null(root.p)) {
    cat("Equal root vector will be used")
    root.p <- rep(1/ncol(Q), ncol(Q))
  }
  
  phy$edge.length[phy$edge.length == 0] <- 1e-05
  
  # F
  
  len.char <- taxa[, 2] %>%
    unique() %>%
    length()
  if (len.char == 1) {
    obj <- NULL
    obj$loglik <- NULL
    obj$diagnostic <- paste("Character is invariant. Analysis stopped.",
                            sep = ""
    )
    return(obj)
  }
  
  
  # workingData <- data.frame(data[, charnum + 1], data[, charnum +
  #                                                       1], row.names = data[, 1])
  # workingData <- workingData[phy$tip.label, ]
  # counts <- table(workingData[, 1])
  # levels <- levels(as.factor(workingData[, 1]))
  # cols <- as.factor(workingData[, 1])
  
  counts <- table(taxa[, 2])
  levels <- names(counts)
  token.maps <- NULL
  
  
  
  cat('\n')
  
  # if(!missing(Edge.Scenarios)){
  #   if(!is.null(Edge.Scenarios$token.maps)){
  #     #token.maps <- Edge.Scenarios$token.maps
  #     token.maps <- Edge.Scenarios$token.maps
  #     
  #     out.tok<-rbind(scen$token.maps, names(scen$token.maps))
  #     out.tok<-data.frame(out.tok)
  #     colnames(out.tok) <- NULL
  #     rownames(out.tok) <-c('original', 'present')
  #     cat('Coding:')
  #     print(out.tok)
  #   }
  # }
  
  cat("State distribution:\n")
  cat("States:", levels, "\n", sep = "\t")
  cat("Counts:", counts, "\n", sep = "\t")
  
  if (!all(colnames(Q) %in%  tk))
    cat(("NOT all states from Q are present in the character.\n"))
  
  if (all(colnames(Q) %in%  tk))
    cat(("ALL states from Q are present in the character.\n"))
  
  
  # factored <- corHMM:::factorData(workingData, charnum = charnum)
  # factored <- My_factor_data(workingData, rate.mat)
  factored <- make_lik_bins(taxa, Qcom = Q)
  nl <- ncol(factored)
  state.names <- colnames(factored)
  # state.names <-colnames(rate.mat)
  bound.hit <- FALSE
  if (ub < 0) {
    ub <- 100
  }
  if (lb < 0) {
    lb <- 0
  }
  if (ub < lb) {
    ub <- 100
    lb <- 0
  }
  
  obj <- NULL
  solution.list <-NULL
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  
  # number of parameters
  np <- lapply(Q, function(x) x[x > 0] %>% unique()) %>%
    unlist() %>%
    unique() %>%
    length()
  #--------------
  liks <- make_lik_bins_all(taxa, Qcom = Q, phy)
  
  lower <- rep(lb, np)
  upper <- rep(ub, np)
  opts <- list(
    algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000",
    ftol_rel = .Machine$double.eps^0.5
  )
  
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### Calculating likelihood using fixed parameters                           ####
  
  if (!is.null(p)) {
    
    cat(
      "\nCalculating likelihood from a set of fixed parameters",
      "\n"
    )
    
    out <- NULL
    # out$solution <- p
    
    out$objective <- likelihood_Sinba(
      p = p, phy = phy, liks = liks, rate = rate, root.p = root.p,
      Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature"
    )
    
    
    loglik <- -out$objective
    est.pars <- out$solution
    
    solution <- matrix(p[rate$Qcom], dim(rate$Qcom))
  }
  
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### Initial Search                                                          ####
  
  #optim.method=c('subplex', 'subplex_annl')
  if (optim.method=='subplex'){
    
  } else if (optim.method=='subplex_annl'){
    
  }
  
  if (is.null(p)) {
    
    
    
    cat('\nStarting', optim.method, 'routine...', "\n")
    cat("\nInitializing...", "\n")
    
    
    if (is.null(ip)) {
      rate.init <- mk_OneRate(Qindv)
      opts <- list(
        algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000",
        ftol_rel = .Machine$double.eps^0.5
      )
      
      dat <- taxa[, c(2)] %>% cbind()
      rownames(dat) <- taxa[, c(1)]
      dat <- phangorn::phyDat(dat, type = "USER", levels = levels(as.factor(dat)))
      
      par.score <- phangorn::parsimony(phy, dat, method = "fitch")
      tl <- sum(phy$edge.length)
      mean.change <- par.score / tl
      
      if (mean.change == 0) {
        ip <- 0.01 + lb
      } else {
        ip <- rexp(1, 1 / mean.change)
      }
      
      if (ip < lb || ip > ub) {
        ip <- lb
      }
    }
    
    
    
    lower.init <- lb
    upper.init <- ub
    
    #-- Subplex
    if (optim.method=='subplex'){
      
      #--- Intials Search
      init <- nloptr(
        x0 = ip,
        eval_f = likelihood_Sinba, lb = lower.init, ub = upper.init,
        opts = opts, phy = phy, liks = liks, rate = rate.init, root.p = root.p,
        Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature"
      )
      
      out <- init
      
      loglik <- -out$objective
      est.pars <- out$solution
      
      solution <- matrix(est.pars[rate.init$Qcom], dim(rate.init$Qcom))
      
      
      ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
      ### Thorough Search                                                         ####
      
      #---- thorough search
      
      lower <- rep(lb, np)
      upper <- rep(ub, np)
      
      if (subplex.rounds==1){
        
        cat(
          "Beginning one thorough search...",
          "\n"
        )
        
        
        out <- nloptr(
          x0 = rep(init$solution, length.out = np), ######## Inference
          eval_f = likelihood_Sinba, lb = lower, ub = upper,
          opts = opts, phy = phy, liks = liks, rate = rate, root.p = root.p,
          Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature"
        )
        
        loglik <- -out$objective
        est.pars <- out$solution
        solution <- matrix(est.pars[rate$Qcom], dim(rate$Qcom))
        
      } else if (subplex.rounds>1){
        
        
        
        cat('Beginning', subplex.rounds, 'thorough searches using', subplex.n.cores, 'cores...\n')
        
        if (mean.change == 0) {
          ip <-runif(subplex.rounds, 0.0001+ lb, 0.02+ lb)
        } else {
          ip <-rexp(subplex.rounds, 1 / mean.change)
        }
        
        if (ip < lb || ip > ub) {
          ip <- rep(lb, subplex.rounds)
        }
        
        
        Lnf <- function(ip){
          out <-nloptr(
            x0 = rep(ip, length.out = np), ######## Inference
            eval_f = likelihood_Sinba, lb = lower, ub = upper,
            opts = opts, phy = phy, liks = liks, rate = rate, root.p = root.p,
            Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature"
          )
          loglik <- -out$objective
          est.pars <- out$solution
          #solution <- matrix(est.pars[rate$Qcom], dim(rate$Qcom))
          list(loglik=loglik, est.pars=est.pars)
        }
        
        out.list <- mclapply(ip, function(x) Lnf(x), mc.cores=subplex.n.cores)
        
        liks.list <- lapply(out.list, function(x) x$loglik) %>% unlist
        liks.best <- which(liks.list==max(liks.list))[1]
        
        solution.list <- lapply(out.list, function(x) x$est.pars) %>% as.data.frame(., col.names =c(1:length(out.list))) %>% rbind(., liks.list)
        loglik <- max(liks.list)
        est.pars <- out.list[[liks.best]]$est.pars
        solution <- matrix(est.pars[rate$Qcom], dim(rate$Qcom))
        
        
      }
      
      #  End thorough search
      # end subplex
      
      
      #--- Subplex annealing
    } else if (optim.method=='subplex_annl'){
      
      lower <- rep(lb, np)
      upper <- rep(ub, np)
      
      cat("Round 1. Beginning simulated annealing...", "\n")
      out.gensa <-GenSA(rep(ip, length.out = np),
                        fn=likelihood_Sinba, lower=lower, upper=upper, control=list(max.call=gensa.its),
                        phy = phy, liks = liks, rate = rate, root.p = root.p, Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature")
      
      
      cat("Round 1. Refining using subplex routine...", "\n")
      out <- nloptr(x0 = out.gensa$par,
                    eval_f = likelihood_Sinba, lb = lower, ub = upper,
                    opts = opts, phy = phy, liks = liks, rate = rate, root.p = root.p,
                    Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature"
      )
      
      
      cat("Round 2. Beginning simulated annealing...", "\n")
      out.gensa <-GenSA(out$solution,
                        fn=likelihood_Sinba, lower=lower, upper=upper, control=list(max.call=gensa.its),
                        phy = phy, liks = liks, rate = rate, root.p = root.p, Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature")
      
      cat("Round 2. Refining using subplex routine...", "\n")
      out <- nloptr(x0 = out.gensa$par,
                    eval_f = likelihood_Sinba, lb = lower, ub = upper,
                    opts = opts, phy = phy, liks = liks, rate = rate, root.p = root.p,
                    Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature"
      )
      
      cat("Round 3. Beginning simulated annealing", "\n")
      out.gensa <-GenSA(out$solution,
                        fn=likelihood_Sinba, lower=lower, upper=upper, control=list(max.call=gensa.its),
                        phy = phy, liks = liks, rate = rate, root.p = root.p, Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature")
      
      cat("Round 3. Refining using subplex", "\n")
      out <- nloptr(x0 = out.gensa$par,
                    eval_f = likelihood_Sinba, lb = lower, ub = upper,
                    opts = opts, phy = phy, liks = liks, rate = rate, root.p = root.p,
                    Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature"
      )
      
      loglik <- -out$objective
      est.pars <- out$solution
      solution <- matrix(est.pars[rate$Qcom], dim(rate$Qcom))
      
    }
    
  }
  
  
  ############################### Diagnostics
  if (diagn == TRUE) {
    
    cat("Finished. Performing diagnostics.", "\n")
    
    
    # stand.error fish.info=H^-1,  se=sqrt(diag(fish.info))
    # conf int [par-(se*1.96), par+(se*1.96)]
    h <- numDeriv::hessian(func = likelihood_Sinba, x = est.pars,
                           phy = phy, liks = liks, rate = rate.init, root.p = root.p,
                           Qindv.zero = Qindv.zero, Pindv = Pindv, edge.pattern.focal = edge.pattern, cub.method = "pcubature")
    
    
    se.vec <- sqrt(diag(corpcor::pseudoinverse(h)))
    up95.lim <- est.pars + (se.vec*1.96)
    low95.lim <- est.pars - (se.vec*1.96)
    solution.se <- matrix(se.vec[rate$Qcom], dim(rate$Qcom))
    solution.up95 <- matrix(up95.lim[rate$Qcom], dim(rate$Qcom))
    solution.low95 <- matrix(low95.lim[rate$Qcom], dim(rate$Qcom))
    rownames(solution.se) <- colnames(solution.se) <- state.names
    rownames(solution.up95) <- colnames(solution.up95) <- state.names
    rownames(solution.low95) <- colnames(solution.low95) <- state.names
    
    hess.eig <- eigen(h, symmetric = TRUE)
    eigval <- signif(hess.eig$values, 2)
    eigvect <- round(hess.eig$vectors, 2)
    
  } else {
    
    solution.se <- NULL
    solution.up95 <- NULL
    solution.low95 <- NULL
    eigval <- NULL
    eigvect <- NULL
  }
  ############################### Making object
  
  if ((any(solution == lb, na.rm = TRUE) || any(solution == ub, na.rm = TRUE)) && (lb != 0 || ub != 100)) {
    bound.hit <- TRUE
  }
  
  rownames(solution) <- colnames(solution) <- state.names
  
  
  obj <- list(
    loglik = loglik, AIC = -2 * loglik + 2 * np,
    AICc = -2 * loglik + (2 * np * (nb.tip / (nb.tip - np - 1))), solution = solution,
    
    opts = opts, data = data,
    phy = phy, iterations = out$iterations,
    bound.hit = bound.hit, n.parameters=np,
    root.p=root.p,
    edge.pattern=edge.pattern,
    Qindv=Qindv, Pindv=Pindv,
    token.maps=token.maps,
    solution.se=solution.se,
    solution.up95=solution.up95,
    solution.low95=solution.low95,
    solution.subplex.list=solution.list,
    optim.method=optim.method,
    eigval= eigval,
    eigvect =eigvect,
    est.pars=est.pars
    
  )
  
  class(obj) <-c(class(obj), 'sinba')
  cat("\\---------------------------------------------------\n")
  
  
  
  return(obj)
}