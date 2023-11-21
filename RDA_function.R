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

doubs_rda$terminfo

ordiplot(doubs_rda, type = "text")
View(rda)
methods(rda)
getAnywhere(rda.default)
getAnywhere(rda.formula)

# Build RDA function ----

x <- vars
y <- species
i=1

y <- as.matrix(y)
envcentre <- colMeans(x)
x <- scale(x, center = envcentre, scale = FALSE)
test2 <- model.matrix(~., as.data.frame(species))[,-1,drop=F]
DISTBASED <- attr(y, "METHOD") == "DISTBASED"
RW <- attr(y, "RW")
CW <- attr(y, "CW")
Q <- qr(x)
?qr
rank <- sum(Q$pivot[seq_len(Q$rank)] > 0)
if (length(Q$pivot) > Q$rank){
  alias <- colnames(Q$qr)[-seq_len(Q$rank)]}else{
  alias <- NULL}
kept <- seq_along(Q$pivot) <= Q$rank & Q$pivot > 0

# So far, we've defined a bunch of null values, and centered and scaled x
# Now, we do QR decomposition
Yfit <- qr.fitted(qr(x), y)
sol <- svd(Yfit)
lambda <- sol$d^2
u <- sol$u
v <- sol$v
## handle zero  eigenvalues and negative eigenvalues
zeroev <- abs(lambda) < max(0, 0 * lambda[1L])
if (any(zeroev)) {
  lambda <- lambda[!zeroev]
  u <- u[, !zeroev, drop = FALSE]}
v <- v[, !zeroev, drop = FALSE]
posev <- lambda > 0
## wa scores
wa <- y %*% v %*% diag(1/sqrt(lambda), sum(posev))


rda_func <- function(x, y){
  # Step 1: Regress species in y over vars in x
  
  preds <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  for (i in 1:ncol(y)) {
    mod <- lm(y[,i] ~ ., data = x)
    preds[,i] <- mod$fitted.values
  }
  colnames(preds) <- colnames(y)
  
  # Step 2: PCA on the fitted values
  pca <- prcomp(preds)
  ordiplot(pca, type = "text")
  biplot(pca)
  pca$sdev
}

## Test function ----
rda_func(x, y)

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

