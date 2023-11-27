#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$     NRES 746 Lecture Script    $#
#$       Intro to Ordination      $#
#$ Martin Genova & Kierstin Acuna $#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# Load libraries ----
library(codep)
library(vegan)
library(datasets)

# Load Doubs fish Data ----

data(Doubs)
species <- as.data.frame(Doubs.fish[-8,])
vars <- as.data.frame(cbind(Doubs.env[-8,],Doubs.geo[-8,]))

?Doubs

# GLM example ----
# Cottus gobio gamma distributed model

CHA.alt <- species$CHA + 1 #get rid of zeros to prep for log link
COGO_mod <- glm(CHA.alt ~ pH + flo + oxy, data = vars, family = Gamma(link = "log"))
summary(COGO_mod) #oxygen is significantly positive!

# Dissimilarity Measures ----

## Types of Distances Coefficients ----

### Euclidean Distance ----

(Y.hmm <- data.frame(hydrophillic_1 = c(1, 0, 0), hydrophillic_2 = c(1, 1, 0),
                    mesic_1 = c(0, 1, 0), mesic_2 = c(0,4,0),
                    xeric_1 = c(0, 1, 3),xeric_2 = c(0, 0, 2), 
                    row.names = c("sample_1_wet", "sample_2_intermediate",
                                                        "sample_3_dry")))

# Calculate Euclidean distance using the dist() function
(Y.hmm.DistEu <- as.matrix(dist(x = Y.hmm, method = "euclidean")))

# Calculate Euclidean Distance by Hand
calc_eu_dist <- function(spe_abun_df) {
  # Create output matrix
  output <- as.data.frame(matrix(NA, nrow = nrow(spe_abun_df), ncol = nrow(spe_abun_df)))
  # Index through the rows of the data frame
  for (i in 1:nrow(spe_abun_df)) {
    x1 <- spe_abun_df[i, ]
    for (t in 1:nrow(spe_abun_df)) {
      x2 <- spe_abun_df[t,]
      # Calculate euclidean distance and place distance into output data frame
      output[i,t] <- sqrt(sum((x1 - x2)^2))
    }
  }
  # Return output
  return(output)
}

# Run Euclidean distance by hand function
(Y.hmm_eu_dist <- calc_eu_dist(Y.hmm))

### Bray-Curtis Coefficient ----

(Y.hmm.BCdist <- vegan::vegdist(Y.hmm, method = "bray", binary = FALSE))

(Y.hmm.BCdist.matrix <- as.matrix(Y.hmm.BCdist))

# Calculate Bray-Curtis coefficient by hand
calc_bc_dist <- function(spe_abun_df) {
  # Create output matrix
  output <- as.data.frame(matrix(NA, nrow = nrow(spe_abun_df), ncol = nrow(spe_abun_df)),
                          row.names = rownames(spe_abun_df))
  colnames(output) <- rownames(spe_abun_df)
  # Index through the rows of the data frame
  for (i in 1:nrow(spe_abun_df)) {
    x1 <- spe_abun_df[i, ]
    for (t in 1:nrow(spe_abun_df)) {
      x2 <- spe_abun_df[t,]
      # Create empty data frame to find the minimum values of each species between two sites
      comp_df <- as.data.frame(matrix(nrow = 2, ncol = ncol(spe_abun_df)))
      # Place the site values into the data frame
      comp_df[1,] = x1
      comp_df[2,] = x2
      # Find the minimum abundance values of each species and sum them.
      min_abundances <- apply(comp_df, 2, min)
      W <- sum(min_abundances)
      # Sum the abundances of site 1
      A = sum(x1)
      # Sum the abundances of site 2
      B = sum(x2)
      # Calculate the Bray-Curtis coefficient
      bc_dist <- (1 - ((2 * W) / (A + B)))
      # Place the BC coefficient into the output data frame
      output[i,t] <- bc_dist
    }
  }
  # Return output
  return(output)
}

# Run Bray-Curtis coefficient by hand function
calc_bc_dist(Y.hmm)


# Unconstrained Ordination ----


## Principal Component Analysis ----


### Breaking Down a PCA ----


# 1.  Load the dataset
data <- Doubs.env[,3:8]

# 2. Standardize the data to have mean zero and unit variance:
(data_std <- as.data.frame(scale(data)))

# 3. Compute the Covariance Matrix
(cov_matrix <- cov(data_std))


# 4. Perform the Eigen-decomposition of the covariance matrix
eigen_decomp <- eigen(cov_matrix)


# Extract Eigenvalues
(eig_values <- eigen_decomp$values)
# Eigenvalues are equal to the sum of squares of the distances of each projected data point in the corresponding principal component
# The sum of squares is maximized in the first principal component.


# Extract Eigenvectors
(eig_vectors <- eigen_decomp$vectors)
rownames(eig_vectors) = colnames(data_std)
colnames(eig_vectors) = c("EV1", "EV2", "EV3", "EV4", "EV5", "EV6")
                          
# Extract the first two eigenvectors
eig_vec_1 <- eig_vectors[,1]
eig_vec_2 <- eig_vectors[,2]


# Note that the first two eigenvectors are perpendicular since the first eigenvector accounts
# for most of the variance in the data and the second eigenvector is orthogonal to the first eigenvector:
eig_vec_1 %*% eig_vec_2


# 5. Show amount of variance contributed by each Principal Component by making a Scree Plot.

  # Calculate the estimated variance for each eigenvalue
(e_var <- eig_values / (nrow(data_std) - 1))

  # Data frame with variance percentages
var_per <- data.frame(
  PC  = c("PC01", "PC02", "PC03", "PC04", "PC05", "PC06"),
  PER = c(e_var) * 100 / sum(e_var) # Calculate the percentage
)
  # Scree plot to show amount of variance accounted by each principal component
barplot(PER ~ PC, data = var_per,
        xlab = "Principal Components",
        ylab = "Percent of Variation %")

# 6. Selecting Principal Components 

  # Kaiser-Guttman Criterion
  # It is generally a good idea to select the principal components that explain most of
  # the variance in the data. One criterion is the Kaiser-Guttman Criterion, which states
  # that any eigenvector with an eigenvalue greater than 1 should be retained.

eig_val_PC <- data.frame(
  PC = c("PC01", "PC02", "PC03", "PC04", "PC05", "PC06"),
  EV = eig_values
)

barplot(EV ~ PC, data = eig_val_PC,
        xlab = "Principal Components",
        ylab = "Eigenvalues"
        )
abline(h = 1, col = "red")

  # Broken stick Model
  # A. Jackson (1993) says that the broken-stick method is one of the better methods for 
  # choosing the number of PCs. The method provides "a good combination of simplicity of 
  # calculation and accurate evaluation of dimensionality relative to the other statistical 
  # approaches" (p. 2212)
  # the "expected proportions" as corresponding to a null model that contains uncorrelated
  # (noise) variables. If you plot the eigenvalues of the correlation matrix against the 
  # broken-stick proportions, the observed proportions that are higher than the expected 
  # proportions indicate which principal components to keep.
  # Fall on conservative side of Principal Component selections.

broken_stick <- function(eig_values) {
  # Calculate Broken Stick Model
  n = length(eig_values)
  bsm = data.frame(j=seq(1:n), prop_var=0)
  bsm$prop_var[1] = 1/n
  for (i in 2:n) {
    bsm$prop_var[i] = bsm$prop_var[i-1] + (1/(n + 1 - i))
    }
  bsm$prop_var = 100*bsm$prop_var/n
  
  # Plot Broken Stick Modol Over 
  barplot(t(cbind(100*eig_values/sum(eig_values), bsm$p[n:1])),
        beside=TRUE, 
        main="Broken Stick Model",
        col=c("red","blue"),
        las=2,
        xlab = "Principal Components", ylab = "Percent of Variation (%)")
  legend("topright", c("Eigenvalues", "Broken stick model"), 
         pch=15,
         col=c("red","blue"), 
         bty="n")

}

broken_stick(eig_values)

# 7. Plot the Principal Components over the data

  # Plot only the first two Principal Components
plot(har ~ pH, col = as.factor(rownames(data_std)), pch = 19,
     xlim=c(-4, 4), ylim = c(-4,4),
     data = (data_std),
     xlab = "pH (Standardized)", ylab = "har (Standardized)")
abline(v=0 , h=0, 
       col = "dark gray")

  # Overlap pertinent eigenvector
abline(0, eig_vec_1[2]/eig_vec_1[1], col='purple')

  # Plot the lines from first eigenvector to points
line1 <- c(0, eig_vec_1[2]/eig_vec_1[1])
perp.segment.coord <- function(x0, y0, line1){
  a <- line1[1]  #intercept 
  b <- line1[2]  #slope
  x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
  y1 <- a + b * x1
  list(x0 = x0, y0 = y0, 
       x1 = x1, y1 = y1)
}
ss <- perp.segment.coord(data_std[,1], data_std[,2], line1)
segments(x0 = ss$x0, x1 = ss$x1, y0 = ss$y0, y1 = ss$y1, col = 'purple')
with(data_std, 
     text(har ~ pH, labels = as.factor(rownames(data_std)), pos = 1, cex=1))

title(main = "First Principal Component over the Standardized Data",
      sub = "Purple Lines Horizontal to the First Principal Components is the Variance", cex.sub = 0.75)

  # Plot both the first and second principal component

  # Make another plot to show second principle component
plot(har ~ pH, col = as.factor(rownames(data_std)), pch = 19,
     xlim=c(-4, 4), ylim = c(-4,4),
     data = (data_std),
     xlab = "pH (Standardized)", ylab = "har (Standardized)")
abline(v=0 , h=0, 
       col = "dark gray")


  #Overlap pertinent eigen-vectors
abline(0, eig_vec_1[2]/eig_vec_1[1], col='purple')
abline(0, eig_vec_2[2]/eig_vec_2[1], col='orange')

# Plot the lines from first eigenvector and second to points
line2 <- c(0, eig_vec_2[2]/eig_vec_2[1])

perp.segment.coord <- function(x0, y0, line2){
  a <- line2[1]  #intercept
  b <- line2[2]  #slope
  x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
  y1 <- a + b * x1
  list(x0 = x0, y0 = y0, 
       x1 = x1, y1 = y1)
}
ss <- perp.segment.coord(data_std[,1], data_std[,2], line2)

segments(x0 = ss$x0, x1 = ss$x1, y0 = ss$y0, y1 = ss$y1,col = 'orange')

with(data_std, text(har ~ pH, labels = as.factor(rownames(data_std)),pos = 1, cex=1))

title(main = "First (Purple) and Second (Orange) Principal Component over the Standardized Data",
      sub = "Lines Horizontal to the Principal Components are the Variance", cex.sub = 0.75)


# 8.Loading Scores

  # Elements of each eigenvector are called loadings and can be interpreted as the contribution of each variable in the data set to 
  # the corresponding principal component.

  # You can make a table with these values and see the contributions of each variable to each principal component:
  
  # Get variable loading scores
variable.loads <- data.frame(
  PC01 = eig_vec_1, # First eigenvector
  PC02 = eig_vec_2  # Second eigenvector
)

  # You can also calculate the loading score for each site, which shows how they are placed in relation to the principal
  # components.

  # Calculate site loading scores
loading.scores <- as.data.frame(as.matrix(data_std) %*% eig_vectors)
colnames(loading.scores) = c("PC01", "PC02", "PC03", "PC04", "PC05", "PC06")

# 9. Make a biplot of the Principal Componants and Variable Loading Scores

  # Set plot parameters
par(mar = c(5, 5, 10, 5),
     mgp = c(2, 1, 0))

plot(loading.scores[,2] ~ loading.scores[,1], 
     xlab = 'PC1', ylab = "PC2",xlim = c(-8,8),ylim = c(-8,8),col = as.factor(rownames(loading.scores)), pch = 19,
     main = "Biplot of Site and Variable Loading Scores Against the First and Second Principal Components")
abline(v = 0, col = "orange")
abline(h = 0, col = "purple")
with(loading.scores, text(PC02 ~ PC01, labels = as.factor(rownames(loading.scores)),pos = 1, cex=1))

par(new=TRUE)

  # Overlay the variable loading scores
plot(PC02 ~ PC01, 
     xlim = c(-0.8, 1), ylim = c(-0.8,1),
     col = "red", pch = 8, axes = F, xlab = "", ylab = "",
     data = variable.loads)
axis(4, ylim = c(-1,1), col = "red")
axis(3, xlim = c(-0.9,1), col = "red")
mtext("Loading Scores",side=3,col="red",line=2.5)  
mtext("Loading Scores",side=4,col="red",line=2.5)  
for (i in 1:nrow(variable.loads)) {
  arrows(x0 = 0, y0 = 0, x1 = variable.loads[i,1],y1 = variable.loads[i,2],
         col = "red", lwd = 1, length = 0.1)
}
with(variable.loads, text(PC02 ~ PC01, labels = as.factor(rownames(variable.loads)),pos = 1, cex=1,
                 col = "red"))
title(main = "Biplot of Site and Variable Loading Scores against the First and Second Principal Components")

### PCA analysis using built-in functions ----

par( mar = c(5, 4, 4, 2) + 0.1,mgp = c(3, 1, 0))

  # Using stats: princomp()
PCA_princomp <- stats::princomp(data_std)
biplot(PCA_princomp)
abline(v= 0, h = 0)
