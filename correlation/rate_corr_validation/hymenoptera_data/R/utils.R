
vectorize_binary <- function(x) {
  c(as.integer(x == 0), as.integer(x == 1))
}


simDirectional_one_char <- function(tree, focal_edge, Q_base, Q_jump){
  focal_edge_length <- tree$edge.length[focal_edge]
  # focal node tipwise
  focal_node <- tree$edge[focal_edge,][2]
  
  #------- Split Tree
  split=list(node=focal_node, bp=0)
  splitTree <- splitTree(tree, split)
  t1 <- splitTree[[1]]
  t2 <- splitTree[[2]]
  #plot(t1)
  #plot(t2)
  
  #------- Sim hist for tree 1
  hist1 <- sim.history(t1, Q_base, nsim=1, direction= "row_to_column")
  #plot(hist1)
  #hist1$tip.label
  #hist1$states
  #hist$tip.label
  
  #------- Sim hist over the focal branch
  init_state <- as.numeric(hist1$states['NA'])
  init_vec <- vectorize_binary(init_state)
  # we need to revert Q if the init_state is 1
  if (init_state==1){
    v <- c(2,1)
    Q_jump <- Q_jump[v,v]
    colnames(Q_jump) <- row.names(Q_jump) <- c(0,1)
    Q_jump
  }
  probs <-init_vec %*% expm::expm(Q_jump*focal_edge_length)
  selected_state <- sample( c(0, 1), size = 1, prob = probs)
  
  #------- Sim hist for tree 2
  hist2 <- sim.history(t2, Q_base, nsim=1, direction= "row_to_column", anc = as.character(selected_state))
  #plot(hist2)
  #hist2$states
  
  #------- Consolidate all tips into one vector
  combined_states <- c(hist1$states, hist2$states)
  tip_states <- combined_states[tree$tip.label]
  # View result
  # hist1$states
  # hist2$states
  # tree$tip.label
  return(tip_states)
}



simDirectional <- function(tree, focal_edge, Q_base, Q_jump, nsim=1){
  # Run simulation N times and collect as list of named vectors
  sims <- replicate(nsim, simDirectional_one_char(tree, focal_edge, Q_base, Q_jump), simplify = FALSE)
  # Combine into a matrix or data frame (tips as rows)
  sim_table <- do.call(cbind, sims)
  return(sim_table)
}
