#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$       NRES 746 Lab Script      $#
#$       Intro to Ordination      $#
#$ Martin Genova & Kierstin Acuna $#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# Load libraries ----
library(codep)
library (vegan)

# Load Doubs fish Data ----
data(Doubs)
species <- as.data.frame(Doubs.fish[-8,])
vars <- as.data.frame(cbind(Doubs.env[-8,],Doubs.geo[-8,]))

# Exercise 1 ----
# Build a function to calculate bray-curtis distance

# Exercise 2 ----
# Compute the likelihood of a model matrix

# Exercise 3 ----
# Build an RDA function

# Answers ----
## Exercise 1 ----
## Exercise 2 ----
## Exercise 3 ----