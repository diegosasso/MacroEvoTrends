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
adj_matrix