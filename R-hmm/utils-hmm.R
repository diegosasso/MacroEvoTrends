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
    
    # Store in environment under character name
    char_name <- paste0("char", charnum)
    assign(char_name, factored, envir = char_env)
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
