# hmm.map <-  c("0", "1")
# hmm.map <-  c("0&1&2", "3&4", "5")
# char_invar_env <-make_invariant_env(data, phy, hmm.map)
# char_invar_env[['char1']]
# char_invar_env[['char2']]
# char_invar_env[['char4']]
make_invariant_env <- function(data, phy, hmm.map) {
  # Merge hmm.map states into one label per character
  hmm.map.merged <- paste0(hmm.map, collapse = "&")
  
  # Build invariant data matrix
  data_invar <- cbind(data[, 1, drop = FALSE], matrix(rep(hmm.map, each = nrow(data)), 
                                                      nrow = nrow(data)))
  
  # Set proper column names: taxa + invariant chars
  colnames(data_invar) <- c(colnames(data)[1], paste0("char_", seq_along(hmm.map)))
  
  # Set the first taxon's invariant chars to the merged mapping
  data_invar[1, 2:ncol(data_invar)] <- hmm.map.merged
  
  # Prepare vectorized invariant characters using your existing function
  char_invar_env <- prepare_vectorized_chars(data_invar, phy, Nchar = length(hmm.map))
  
  # Update the first row of each invariant character likelihood matrix
  char_names <- ls(char_invar_env)
  for (ci in seq_along(char_names)) {
    char <- char_names[ci]
    char_invar_env[[char]][1, ] <- char_invar_env[[char]][2, ]
  }
  
  return(char_invar_env)
}



prepare_vectorized_chars <- function(data, phy, Nchar) {
  # Create an environment for ultrafast constant-time access
  char_env <- new.env(hash = TRUE, parent = emptyenv())
  
  # Loop over characters
  for (charnum in 1:Nchar) {
    # Prepare input for factorData
    data.sort <- data.frame(data[, charnum + 1], data[, charnum + 1],
                            row.names = data[, 1])
    data.sort <- data.sort[phy$tip.label, ]
    
    # Factorize using corHMM
    factored <- corHMM:::factorData(data.sort, charnum = charnum)
    model.set.final <- corHMM:::rate.cat.set.rayDISC(phy = phy, data = data.sort, model = 'ER', charnum = 1)
    
    # Store in environment under character name
    char_name <- paste0("char", charnum)
    assign(char_name, model.set.final$lik, envir = char_env)
  }
  
  return(char_env)
}

char_env_summary <- function(char_env) {
  # Iterate through all characters in the environment
  for (char_name in ls(char_env)) {
    mat <- char_env[[char_name]]
    
    # Number of taxa = number of rows
    n_taxa <- nrow(mat)
    
    # Number of possible states = number of columns
    n_states <- ncol(mat)
    
    # Count how many taxa are in each state
    # Pick the column names if they exist, otherwise use numeric labels
    state_labels <- if (!is.null(colnames(mat))) colnames(mat) else seq_len(n_states)
    
    # Find the max column (one-hot encoding = 1 for observed state)
    state_counts <- colSums(mat)
    names(state_counts) <- state_labels
    
    # Print summary
    cat("Character:", char_name, "\n")
    cat("  Number of taxa:  ", n_taxa, "\n")
    cat("  Number of states:", n_states, "\n")
    cat("  State counts:\n")
    print(state_counts)
    cat("-------------------------------\n")
  }
}
