library(codep)
library(DHARMa)

# Load Doubs fish Data ----

data(Doubs)
species <- as.data.frame(Doubs.fish[-8,])
vars <- as.data.frame(cbind(Doubs.env[-8,],Doubs.geo[-8,]))

# Simple linear regression
CHA.alt <- species$CHA + 1
COGO_mod <- glm(CHA.alt ~ pH + flo + oxy, data = vars, family = Gamma(link = "log"))
summary(COGO_mod)
