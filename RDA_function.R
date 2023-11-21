# Set working directory ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list.files()

# Load Libraries ----

library(codep)
library(vegan)

# Load Doubs fish Data ----

data(Doubs)
species <- as.data.frame(Doubs.fish[-8,])
vars <- as.data.frame(cbind(Doubs.env[-8,],Doubs.geo[-8,]))

# Run black box RDA ----

# To simplify things, pick 3 variables to work with:
#   flo, oxy, and pH

vars <- vars[,c("flo", "oxy", "pH")]

species <- decostand(species, method = "hellinger")
# The hellinger transformation helps to account for the double zero problem, it's the 
#   square root of the relative abundance at each site.
# If I have them do a PCA, the data should be hellinger-transformed. But if I have them
#   do a PCoA, they can leave the data as is, because the bray-curtis distance metric also
#   accounts for the double-zero problem.

vars <- decostand(vars, method = "standardize")
# Continuous variables need to be scaled and centered

doubs_rda <- rda(species ~ ., data = vars)

summary(doubs_rda)
# 47% of the variance is constrained, meaning 47% of the variance in Y is explained
#   by the environmental variables in X

ordiplot(doubs_rda, type = "text")


# Build RDA function ----

x <- vars
y <- species
i=1

rda_func <- function(x, y){
  # Step 1: Regress species in y over vars in x
  
  preds <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  for (i in 1:ncol(y)) {
    mod <- lm(y[,i] ~ ., data = x)
    preds[,i] <- mod$fitted.values
  }
  colnames(preds) <- colnames(y)
 summary(mod) 
  
  # Step 2: PCA on the fitted values
  pca <- prcomp(preds)
  biplot(pca)
  
}

## Test function ----


# Old function building on dune data
# Inputs: x and y matrices

rda_func <- function(x, y){
  # Step 1: Regress species in y over vars in x
  preds <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  for (i in 1:ncol(y)) {
    mod <- lm(y[,i] ~ x[,1] + x[,2] + x[,3] + x[,4] + x[,5])
    preds[,i] <- mod$fitted.values
  }
  
  # Step 2: PCA on the fitted values
  pca <- prcomp(preds)
  biplot(pca)
  
}

