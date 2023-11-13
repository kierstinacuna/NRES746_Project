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

# Explore data ----

hist(species$CHA)

# Run black box RDA ----

species <- decostand(species, method = "hellinger")

vars <- decostand(vars, method = "standardize")

doubs_rda <- rda(species ~ ., data = vars)

summary(doubs_rda)

prcomp(doubs_rda)

ordiplot(doubs_rda,
         type = "text"
         )


# Build RDA function ----

rda_func <- function(x, y){
  # Step 1: Regress species in y over vars in x
  preds <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  for (i in 1:ncol(y)) {
    mod <- lm(y[,i] ~ x[,1] + x[,2] + x[,3] + x[,4] + x[,5] + x[,6] + x[,7] + x[,8] + x[,9] +
                x[,10] + x[,11])
    preds[,i] <- mod$fitted.values
  }
  
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

