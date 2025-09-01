# compute_aic <- function(file_path, ...) {
#   # Load dataset
#   data_hmm <- readRDS(file_path)
#   
#   # Extract useful info
#   Ndata <- length(data_hmm)
#   Nchar <- ncol(data_hmm[[1]])-1
#   dataset_name <- basename(file_path)
#   
#   # Prepare storage for AIC
#   AIC_mat <- matrix(NA, nrow = Ndata, ncol = 3)
#   colnames(AIC_mat) <- c("Q.r1", "Q.r2.asym", "Q.r2.sym")
#   rownames(AIC_mat) <- paste0(dataset_name, "_sim", seq_len(Ndata))
#   
#   # Loop over all simulations in the file
#   for (i in seq_len(Ndata)) {
#     fit_Q.r1 <- rayDISC_multi(tree, data_hmm[[i]], Nchar = Nchar,
#                               rate.mat = Q.r1, hmm.map = c("0&1", "2&3"),
#                               root.p = "flat", node.states = "none",
#                               lewis.asc.bias = TRUE)
#     
#     fit_Q.r2.asym <- rayDISC_multi(tree, data_hmm[[i]], Nchar = Nchar,
#                                    rate.mat = Q.r2.asym, hmm.map = c("0&1", "2&3"),
#                                    root.p = "flat", node.states = "none",
#                                    lewis.asc.bias = TRUE)
#     
#     fit_Q.r2.sym <- rayDISC_multi(tree, data_hmm[[i]], Nchar = Nchar,
#                                   rate.mat = Q.r2.sym, hmm.map = c("0&1", "2&3"),
#                                   root.p = "flat", node.states = "none",
#                                   lewis.asc.bias = TRUE)
#     
#     # Save AIC values
#     AIC_mat[i, ] <- c(fit_Q.r1$AIC, fit_Q.r2.asym$AIC, fit_Q.r2.sym$AIC)
#   }
#   
#   return(AIC_mat)
# }


compute_aic <- function(file_path) {
  # Load dataset
  data_hmm <- readRDS(file_path)
  
  # Extract useful info
  Ndata <- length(data_hmm)
  Nchar <- ncol(data_hmm[[1]])-1
  dataset_name <- basename(file_path)
  
  # Start a local cluster
  ncores <- max(1, detectCores() - 1)
  cl <- makeCluster(ncores)
  
  # Source helper files on each worker
  clusterEvalQ(cl, {
    source("R-hmm/utils-hmm.R")
    source("R-hmm/rayDISC-multi.R")
  })
  
  # Export needed objects to workers
  clusterExport(
    cl,
    varlist = c("tree", "data_hmm", "Nchar", 
                "Q.r1", "Q.r2.asym", "Q.r2.sym", "rayDISC_multi"),
    envir = environment()
  )
  
  # Prepare storage for AIC
  # AIC_mat <- matrix(NA, nrow = Ndata, ncol = 3)
  # colnames(AIC_mat) <- c("Q.r1", "Q.r2.asym", "Q.r2.sym")
  # rownames(AIC_mat) <- paste0(dataset_name, "_sim", seq_len(Ndata))
  
  #for (i in seq_len(Ndata)) {
  # Parallelize over simulations in this file
  AIC_list <- parLapply(cl, seq_len(Ndata), function(i) {
    fit_Q.r1 <- rayDISC_multi(tree, data_hmm[[i]], Nchar = Nchar,
                              rate.mat = Q.r1, hmm.map = c("0&1", "2&3"),
                              root.p = "flat", node.states = "none",
                              lewis.asc.bias = TRUE)
    
    fit_Q.r2.asym <- rayDISC_multi(tree, data_hmm[[i]], Nchar = Nchar,
                                   rate.mat = Q.r2.asym, hmm.map = c("0&1", "2&3"),
                                   root.p = "flat", node.states = "none",
                                   lewis.asc.bias = TRUE)
    
    fit_Q.r2.sym <- rayDISC_multi(tree, data_hmm[[i]], Nchar = Nchar,
                                  rate.mat = Q.r2.sym, hmm.map = c("0&1", "2&3"),
                                  root.p = "flat", node.states = "none",
                                  lewis.asc.bias = TRUE)
    
    # Save AIC values
    #AIC_mat[i, ] <- c(fit_Q.r1$AIC, fit_Q.r2.asym$AIC, fit_Q.r2.sym$AIC)
    return(c(fit_Q.r1$AIC, fit_Q.r2.asym$AIC, fit_Q.r2.sym$AIC))
  })
  
  AIC_mat <- do.call(rbind, AIC_list)
  colnames(AIC_mat) <- c("Q.r1", "Q.r2.asym", "Q.r2.sym")
  rownames(AIC_mat) <- paste0(dataset_name, "_sim", seq_len(Ndata))
  
  return(AIC_mat)
}
