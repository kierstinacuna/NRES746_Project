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

## Explore the data ----
# Count the number of species frequencies in each abundance class
ab <- table(unlist(species))
# Plot distribution of species frequencies
barplot(ab, xlab = "Abundance class", ylab = "Frequency", col = grey(5:0/5))

# Look at the spread of the predictor variables
?Doubs
summary(vars)

# Exercise 1 ----
# Build a function to prepare the data for RDA
RDA_prep <- function(x, y){
  # Step 1: hellinger-transform the y data
  
  # Step 2: center and scale the x data
}

# Exercise 2 ----
# Compute the likelihood of a model matrix

# Exercise 3 ----
# Build an RDA function

# Answers ----
## Exercise 1 ----
x <- vars
y <- species
i <- 1

RDA_prep <- function(x, y){
 
  # Step 1: hellinger-transform the y data
  yhell <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  y <- as.matrix(y)
  for (i in 1:nrow(y)) {
    yhell[i,] <- sqrt(y[i,]/sum(y[i,]))
  }
  colnames(yhell) <- colnames(y)

  # Step 2: center and scale the x data
  xscale <- scale(x)
  
  # Output
  out <- list()
  out[["x"]] <- as.data.frame(xscale)
  out[["y"]] <- as.data.frame(yhell) 
  return(out)
}

data <- RDA_prep(vars, species)
data$x
data$y

## Exercise 2 ----
## Exercise 3 ----