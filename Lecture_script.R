####################################
##     NRES 746 Lecture Script    ##
##       Intro to Ordination      ##
## Martin Genova & Kierstin Acuna ##
####################################

# Load libraries ----
library(codep)

# Load Doubs fish Data ----

data(Doubs)
species <- as.data.frame(Doubs.fish[-8,])
vars <- as.data.frame(cbind(Doubs.env[-8,],Doubs.geo[-8,]))

# GLM example ----
# Cottus gobio gamma distributed model

CHA.alt <- species$CHA + 1 #get rid of zeros to prep for log link
COGO_mod <- glm(CHA.alt ~ pH + flo + oxy, data = vars, family = Gamma(link = "log"))
summary(COGO_mod) #oxygen is significantly positive!
