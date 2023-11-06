# Final Project Testing Script

# Packages

library(tidyverse)
library(vegan)

# Possible Data Sets

# Dune vegetation and environment data from Dutch Dune Meadows
?dune
data(dune)
data(dune.env)

hist(dune$Agrostol, freq = F)

# JAGS File
cat(
  "model {
  
  # Likelihood
  
  for (i in 1:n.obs) {
  exp_cover[i] <- Use[Management[Manure[Moisture[i]]]] + b.a1 * A1[i]
  
  p_cover.class[i] ~ dnegbin(p = , r = )
  
  }
  
  # Prior
  "
)
