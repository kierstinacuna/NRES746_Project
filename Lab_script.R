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

# Use the function to prepare our data for RDA
data <- RDA_prep(         )

# Exercise 2 ----
# Step 1: Create a global RDA using all the predictor variables, and a null RDA 
#   using none of the predictor variables.
rda_global <- rda(          )

rda_null <- rda(          )  

# Step 2: Use the ordiR2step function for variable selection. Define the object as the 
#   null rda, and the scope as the global rda.
selection <- ordiR2step(         )

# Step 3: View the results. The formula in the "Call" contains the variables that 
#   ordiR2step selected.
selection

# Exercise 3 ----


# Answers ----
## Exercise 1 ----
RDA_prep <- function(x, y){
 
  # Step 1: hellinger-transform the y data
  yhell <- matrix(NA, nrow = nrow(y), ncol = ncol(y)) #create a storage matrix
  y <- as.matrix(y) #convert y to a matrix
  for (i in 1:nrow(y)) {
    yhell[i,] <- sqrt(y[i,]/sum(y[i,])) #use the hellinger-transformation function
  }
  colnames(yhell) <- colnames(y) #assign the species names as column names to the new
                                 #  hellinger-transformed matrix

  # Step 2: center and scale the x data
  xscale <- scale(x)
  
  # Output
  out <- list()
  out[["x"]] <- as.data.frame(xscale)
  out[["y"]] <- as.data.frame(yhell) 
  return(out)
}

# Compare our answer to the vegan::decostand answer
our_hell <- RDA_prep(vars, species)$y
vegan_hell <- decostand(species, "hellinger")
which(our_hell != vegan_hell) #they're the same!

# Use our function to prepare our data for RDA
data <- RDA_prep(vars, species)

## Exercise 2 ----
# Step 1: Use the RDA_prep function to prepare the data. Then, create a global RDA 
#   using all the predictor variables, and a null RDA using none of the predictor variables.

rda_global <- rda(data$y ~ ., data=data$x)

rda_null <- rda(data$y ~ 1, data=data$x)

# Step 2: Use the ordiR2step function for variable selection
selection <- ordiR2step(object = rda_null, scope = rda_global)

# Step 3: View the results. The formula in the "Call" contains the variables that 
#   ordiR2step selected.
selection

## Exercise 3 ----

