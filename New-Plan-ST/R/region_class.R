# Body regions
regions <- c(
  "cranium", "mouthparts", "pronotum", "propectus", 
  "mesonotum", "mesopectus", "metanotum", "metapectal-propodeal complex",
  "metasoma", "female genitalia", 
  "fore leg", "mid leg", "hind leg", 
  "fore wing", "hind wing", "sclerites",      "muscles"
)

# Classification
groups <- c(
  "head", "head",          # cranium, mouthparts
  "mesosoma", "mesosoma",  # pronotum, propectus
  "mesosoma", "mesosoma", "mesosoma", "mesosoma",  # mesonotum...complex
  "metasoma", "metasoma",  # metasoma, female genitalia
  "legs", "legs", "legs",          # legs
  "wings", "wings",                 # wings
  "sclerites",  "muscles"
)



# Create object
region_class <- data.frame(
  region = regions,
  group = groups,
  stringsAsFactors = FALSE
)

#region_class
