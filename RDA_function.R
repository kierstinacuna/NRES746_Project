# Set working directory ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list.files()

# Download Data ----

## Dune Data ----
library(vegan)
data(dune)
data(dune.env)

## Doubs Data ----
species <- read.csv("workshop10-main/book-en/data/doubsspe.csv", row.names = 1)
species <- species[-8,]

vars <- read.csv("workshop10-main/book-en/data/doubsenv.csv", row.names = 1)
vars <- vars[-8,]

# Explore data ----

## Dune data ----
hist(dune$Achimill)
hist(dune$Agrostol)
hist(dune$Airaprae)
hist(dune$Alopgeni)

## Doubs data ----
hist(species$CHA)


# Run black box RDA ----

# Dune
ordiplot(rda(dune ~ ., data=dune.env)
         #, type = "text"
         )

# Doubs

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



# Inputs: x and y matrices

x <- dune.env
y <- dune
i <- 1

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

rda_func(x = dune.env, y = dune)

rda(dune ~ ., dune.env)

anova(mod)
rda
pca$x
pca$rotation
pca$x
summary(pca)
summary(mod)

