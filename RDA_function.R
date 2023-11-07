library(vegan)
data(dune)
data(dune.env)

hist(dune$Achimill)
hist(dune$Agrostol)
hist(dune$Airaprae)
hist(dune$Alopgeni)

ordiplot(rda(dune ~ ., data=dune.env), type = "text")

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
anova(mod)
rda
pca$x
pca$rotation
summary(pca)
summary(mod)
