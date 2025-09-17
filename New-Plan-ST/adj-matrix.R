# Define body parts
parts <- c(
  "cranium", "mouthparts", "pronotum", "propectus", "mesonotum", "mesopectus",
  "metanotum", "metapectal-propodeal complex", "metasoma", "female genitalia",
  "fore leg", "mid leg", "hind leg", "fore wing", "hind wing"
)

# Initialize adjacency matrix
adj_matrix <- matrix(0, nrow = length(parts), ncol = length(parts))
rownames(adj_matrix) <- parts
colnames(adj_matrix) <- parts

# Define pairwise adjacencies based on topological contact
adj_pairs <- list(
  c("cranium", "mouthparts"),
  c("cranium", "pronotum"),
  c("pronotum", "propectus"),
  c("propectus", "mesopectus"),
  c("mesopectus", "mesonotum"),
  c("mesonotum", "metanotum"),
  c("metanotum", "metapectal-propodeal complex"),
  c("metapectal-propodeal complex", "metasoma"),
  c("metasoma", "female genitalia"),
  c("propectus", "fore leg"),
  c("mesopectus", "mid leg"),
  c("metapectal-propodeal complex", "hind leg"),
  c("mesonotum", "fore wing"),
  c("metanotum", "hind wing")
)

# Fill in the adjacency matrix (symmetric)
for (pair in adj_pairs) {
  i <- pair[1]
  j <- pair[2]
  adj_matrix[i, j] <- 1
  adj_matrix[j, i] <- 1
}

# View matrix
#adj_matrix

#---------------



# 1. Assign body parts to anatomical groups
group_map <- list(
  head         = c("cranium", "mouthparts"),
  mesosoma     = c("pronotum", "propectus", "mesonotum", "mesopectus", "metanotum", "metapectal-propodeal complex"),
  metasoma_gen = c("metasoma", "female genitalia"),
  wings        = c("fore wing", "hind wing"),
  fore_leg     = "fore leg",
  mid_leg      = "mid leg",
  hind_leg     = "hind leg"
)

# Flatten group map into a named vector: body part → group
body_parts <- unlist(group_map)
part_to_group <- rep(names(group_map), lengths(group_map))
names(part_to_group) <- body_parts

# 2. Define adjacency between major anatomical groups
group_edges <- data.frame(
  from = c("head",     "mesosoma",      "mesosoma",      "mesosoma",   "mesosoma",   "mesosoma"),
  to   = c("mesosoma", "metasoma_gen",  "fore_leg",      "mid_leg",    "hind_leg",   "wings")
)

# Build group graph
group_graph <- graph_from_data_frame(group_edges, directed = FALSE)
group_distances <- distances(group_graph)

# 3. All pairwise combinations of body parts
all_parts <- names(part_to_group)
all_pairs <- combn(all_parts, 2, simplify = FALSE)

# 4. Assign distances
get_distance <- function(p1, p2) {
  g1 <- part_to_group[p1]
  g2 <- part_to_group[p2]
  if (g1 == g2) {
    return(1)
  } else {
    return(1 + group_distances[g1, g2])  # 1 + number of inter-group steps
  }
}

# Compute distances
edges <- lapply(all_pairs, function(pair) {
  dist <- get_distance(pair[1], pair[2])
  c(pair[1], pair[2], dist)
})

# Format as data frame
edge_df <- as.data.frame(do.call(rbind, edges))
colnames(edge_df) <- c("from", "to", "weight")
edge_df$weight <- as.numeric(edge_df$weight)

# 5. Build graph from all body parts
g2 <- graph_from_data_frame(edge_df, directed = FALSE)

# Optional: visualize
# plot(g,
#      edge.label = E(g)$weight,
#      layout = layout_with_fr(g),
#      vertex.label.cex = 0.8,
#      vertex.size = 30,
#      edge.width = 1,
#      edge.color = "gray")

# 6. Get distance matrix (optional)
#body_part_distance_matrix <- distances(g, weights = E(g)$weight)



#-------------

ana.d <- c(
  "cranium vs pronotum" = 2,
  "cranium vs propectus" = 2,
  "cranium vs mesonotum" = 2,
  "cranium vs metapectal-propodeal complex" = 2,
  "mouthparts vs mesonotum" = 2,
  "mouthparts vs mesopectus" = 2,
  "mouthparts vs metapectal-propodeal complex" = 2,
  "pronotum vs propectus" = 1,
  "pronotum vs mesonotum" = 1,
  "pronotum vs mesopectus" = 1,
  "pronotum vs metapectal-propodeal complex" = 1,
  "propectus vs mesonotum" = 1,
  "propectus vs mesopectus" = 1,
  "propectus vs metapectal-propodeal complex" = 1,
  "mesonotum vs mesopectus" = 1,
  "mesonotum vs metapectal-propodeal complex" = 1,
  "mesopectus vs metapectal-propodeal complex" = 1,
  "mesopectus vs metasoma" = 2,
  "metanotum vs metasoma" = 2,
  "propectus vs metasoma" = 2,
  "mesopectus vs hind wing" = 2,
  "mesonotum vs fore wing" = 2,
  "mesonotum vs hind wing" = 2,
  "fore wing vs hind wing" = 1,
  "metapectal-propodeal complex vs mid leg" = 2,
  "pronotum vs mid leg" = 1
)


# Body regions
regions <- c(
  "cranium", "mouthparts", "pronotum", "propectus", 
  "mesonotum", "mesopectus", "metanotum", "metapectal-propodeal complex",
  "metasoma", "female genitalia", 
  "fore leg", "mid leg", "hind leg", 
  "fore wing", "hind wing"
)

# Classification
groups <- c(
  "head", "head",          # cranium, mouthparts
  "mesosoma", "mesosoma",  # pronotum, propectus
  "mesosoma", "mesosoma", "mesosoma", "mesosoma",  # mesonotum...complex
  "metasoma_gen", "metasoma_gen",  # metasoma, female genitalia
  "legs", "legs", "legs",          # legs
  "wings", "wings"                 # wings
)

# Create object
region_class <- data.frame(
  region = regions,
  group = groups,
  stringsAsFactors = FALSE
)

region_class
