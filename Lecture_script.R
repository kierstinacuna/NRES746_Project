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

# GLM example ----
# Cottus gobio gamma distributed model

CHA.alt <- species$CHA + 1 #get rid of zeros to prep for log link
COGO_mod <- glm(CHA.alt ~ pH + flo + oxy, data = vars, family = Gamma(link = "log"))
summary(COGO_mod) #oxygen is significantly positive!

# Dissimilarity Measures ----

## Types of Distances Coefficients ----

### Euclidean Distance ----

Y.hmm <- data.frame(hydrophillic_1 = c(1, 0, 0), hydrophillic_2 = c(1, 1, 0),
                    mesic_1 = c(0, 1, 0), mesic_2 = c(0,4,0),
                    xeric_1 = c(0, 1, 3),xeric_2 = c(0, 0, 2), 
                    row.names = c("sample_1_wet", "sample_2_intermediate",
                                                        "sample_3_dry"))

# Calculate Euclidean distance using the dist() function
(Y.hmm.DistEu <- dist(x = Y.hmm, method = "euclidean"))

# Coerce dist data type into a matrix
(eu_dist_matrix_y.hmm <- as.matrix(Y.hmm.DistEu))

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

calc_bc_dist(Y.hmm)


# Unconstrained Ordination ----

## Principal Component Analysis ----

### Breaking Down a PCA ----

# 1.  Load the varechem dataset
data(varechem)
(data <- varechem[, 1:2]) # Select the Nitrogen and Phosphous data columns

# 2. Standardize the data to have mean zero and unit variance:
data_std <- scale(data)

# 3. Compute the Covariance Matrix
(cov_matrix <- cov(data_std))

# 4. Perform the Eigen-decomposition of the covariance matrix
eigen_decomp <- eigen(cov_matrix)
Eigenvalues <- eigen_decomp$values
Eigenvectors <- eigen_decomp$vectors

# Function to plot the eigenvectors in relation to the data
plot_2pc <- function(eig_vec, eig_val, std_data) {
  # Set Up Data and empty vectors
  eig_vec. <- as.data.frame(eig_vec)
  row.names(eig_vec.) <- c("P", "N")
  colnames(eig_vec.) <- c("PC1", "PC2")
  eig_vec.[, 1] <- eig_vec.[, 1]*-1
  Y_std <- as.data.frame(std_data)
  
  # Set Plotting Parameters
  op <- par(mfrow = c(2, 1),     # 2x2 layout
            oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
            mar = c(1, 1, 0, 0), # space for one row of text at ticks and to separate plots
            mgp = c(2, 1, 0)    # axis label at 2 rows distance, tick labels at 1 row
  )       
  
  # Plot P against N to show first principal component
  plot(N ~ P, col = as.factor(rownames(Y_std)), pch = 19, xlim=c(-2.2, 2.2), 
       ylim = c(-2.2,2.2), data = as.data.frame(Y_std))
  abline(v=0 , h=0, 
         col = "dark gray")
  
  # Overlap pertinent eigenvector
  abline(0, eig_vec[2, 1]/eig_vec[1, 1], col='purple')
  arrows(x0 = 0, y0 = 0, x1 = eig_val[1]*eig_vec[1, 1],y1 = eig_val[1]*eig_vec[2, 1],
         col = "purple", lwd = 2)
  
  # Plot the lines from first eigenvector to points
  line1 <- c(0, eig_vec[2, 1]/eig_vec[1, 1])
  perp.segment.coord <- function(x0, y0, line1){
    a <- line1[1]  #intercept 
    b <- line1[2]  #slope
    x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
    y1 <- a + b * x1
    list(x0 = x0, y0 = y0, 
         x1 = x1, y1 = y1)
  }
  ss <- perp.segment.coord(Y_std$P, Y_std$N, line1)
  
  # do.call(segments, ss)
  # which is the same as:
  segments(x0 = ss$x0, x1 = ss$x1, y0 = ss$y0, y1 = ss$y1, col = 'purple')
  points(N ~ P, col = as.factor(rownames(Y_std)), pch = 19, data = Y_std)
  with(Y_std, 
       text(N ~ P, labels = as.factor(rownames(Y_std)), pos = 1, cex=1.4))
  
  # Make another plot to show second principle component
  plot(N ~ P, col = as.factor(rownames(Y_std)), pch = 19, xlim=c(-2.2, 2.2), 
       ylim = c(-2.2,2.2), data = as.data.frame(Y_std))
  abline(v=0 , h=0, col = "dark gray")
  
  #Overlap pertinent eigen-vectors
  abline(0, eig_vec[2, 1]/eig_vec[1, 1], col='purple')
  abline(0, eig_vec[1, 2]/eig_vec[2, 2], col='orange')
  
  arrows(x0 = 0, y0 = 0, x1 = eig_val[1]*eig_vec[1, 1], y1 = eig_val[1]*eig_vec[2, 1], col = "purple", lwd = 2)
  
  arrows(x0 = 0,y0 = 0, x1 = eig_val[2]*eig_vec[1,2], y1 = eig_val[2]*eig_vec[2, 2], col = "orange", lwd = 2)
  
  line2 <- c(0, eig_vec[1, 2]/eig_vec[1, 1])
  
  perp.segment.coord <- function(x0, y0, line2){
    a <- line2[1]  #intercept
    b <- line2[2]  #slope
    x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
    y1 <- a + b * x1
    list(x0 = x0, y0 = y0, 
         x1 = x1, y1 = y1)
  }
  ss <- perp.segment.coord(Y_std$P, Y_std$N, line2)
  
  segments(x0 = ss$x0, x1 = ss$x1, y0 = ss$y0, y1 = ss$y1,col = 'orange')
  
  points(N ~ P, col = as.factor(rownames(Y_std)), pch = 19, data = Y_std)
  
  with(Y_std, text(N ~ P, labels = as.factor(rownames(Y_std)),pos = 1, cex=1.4))
  
  title(xlab = "N",ylab = "P", outer = TRUE, line = 3)
}

plot_2pc(Eigenvectors, Eigenvalues, data_std)

# 5. Project the standardized data onto the Eigen-space
F_PrComps <- data_std %*% Eigenvectors
head(F_PrComps)


# Trying to plot the eigen vectors


plot(data_std[,1] ~ data_std[,2])

ggplot(data.frame(data_set_2), aes(x = varc_1, y = varc_2)) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = 0, size = .5) +
  geom_hline(yintercept = 0, size = .5) +
  geom_abline(slope = ev1_m, color = "blue", size = 0.7) +
  geom_abline(slope = ev2_m, color = "red", size = 0.7) +
  theme_classic()


### Using blackbox PCA functions ----

data(varechem)

Y <- varechem[, 1:2]
Y_std <- as.matrix(scale(Y))
Y_R <- cov(Y_std)

(Eigenvalues <- eigen(Y_R, symmetric = T)$values)
(Eigenvectors <- eigen(Y_R, symmetric = T)$vectors)

(F_PrComps <- Y_std %*% Eigenvectors)

head(F_PrComps)

# Using stats::prcom()
PCA_prcomp <- stats::prcomp(Y, center = TRUE, scale = TRUE)

# or PCA_prcomp <- prcomp(Y_std)

(prcomp_scores <- PCA_prcomp$x)

# Using stats: princomp()

PCA_princomp <- stats::princomp(Y_std)

head(PCA_princomp$scores)

# Using vegan::rda()

PCA_vegan_rda <- vegan::rda(Y_std)

test <- scores(PCA_vegan_rda, display = "sites", scaling = 1, choices = seq_len(PCA_vegan_rda$CA$rank),
       const = sqrt(PCA_vegan_rda$tot.chi * (nrow(PCA_vegan_rda$CA$u) -
                                               1)))[1:5, ]


### Practice using on Hellinger transformed Doubs fish data ----

# Hellinger transform Doubs fish data
spe.hel <- vegan::decostand(species, method = "hellinger")

# Standardize Hellinger transformed fish data
spe.hel.std <- as.matrix(scale(spe.hel))

# Calculate covariance matrix of the fish data
spe.cov <- cov(spe.hel.std)

# Calculate the eigenvalues and eigenvectors
spe_cov_e <- eigen(spe.cov)
(spe.eigvec <- spe_cov_e$vectors)
(spe.eigval <- spe_cov_e$values)


(spe.pcscores <- spe.hel.std %*% spe.eigvec)

spe.pca.vegan <- rda(spe.hel)
summary(spe.pca.vegan)
plot(spe.pca.vegan$CA$eig, xlim = c(0,50))
points(spe.eigval, col = "red")
  


# From R-bloggers:
# Calculate PCA from scratch

library(tidyverse)

set.seed(1)
# Variable 1
var_1 <- rnorm(50, 50, sd = 3)
# Variable 2
var_2 <- .5*var_1 + rnorm(50, sd = sqrt(3))
# Both variables in a data frame
data_set_1 <- data.frame(var_1, var_2)
head(data_set_1)

# A scatter plot with the two simulated variables
ggplot(data_set_1, aes(x = var_1, y = var_2)) +
  geom_point(color = "blue", size = 2) +
  xlab("Variable 1") +
  ylab("Variable 2") +
  theme_classic()

# Subtract the mean from each variable
data_set_1 <- data_set_1 %>% 
  mutate(varc_1 = var_1 - mean(var_1), varc_2 = var_2 - mean(var_2))
head(data_set_1)

# Scatter plot for the centered data 
ggplot(data_set_1, aes(x = varc_1, y = varc_2)) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = 0, size = .5) +
  geom_hline(yintercept = 0, size = .5) +
  theme_classic()

# Select just the centered variables
data_set_2 <- data_set_1 %>% 
  select(varc_1, varc_2) %>% 
  as.matrix()
# Calculate the covariance matrix
cov_m <- (t(data_set_2) %*% data_set_2) / (nrow(data_set_2) - 1) 
cov_m

# Use eigen() to obtain eigenvectors and eigenvalues
cov_e <- eigen(cov_m)
# Eigenvectors
e_vec <- cov_e$vectors
# Eigenvalues
e_val <- cov_e$values

# First eigenvector 
ev_1 <- e_vec[,1]
# First eigenvector slope
ev1_m <- ev_1[2] / ev_1[1]
# Second eigenvector 
ev_2 <- e_vec[,2]
# Second eigenvector slope
ev2_m <- ev_2[2] / ev_2[1]
# Scatter plot showing the span of both eigenvectors 
ggplot(data.frame(data_set_2), aes(x = varc_1, y = varc_2)) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = 0, size = .5) +
  geom_hline(yintercept = 0, size = .5) +
  geom_abline(slope = ev1_m, color = "blue", size = 0.7) +
  geom_abline(slope = ev2_m, color = "red", size = 0.7) +
  theme_classic()

# Multiply both eigenvectors 
ev_1 %*% ev_2
##      [,1]
## [1,]    0

# Calculate the estimated variance for each eigenvalue
e_var <- e_val / (nrow(data_set_2) - 1)
# Data frame with variance percentages
var_per <- data.frame(
  PC  = c("PC1", "PC2"),
  PER = c(e_var) * 100 / sum(e_var) # Calculate the percentage
)
# Scree plot 
ggplot(var_per, aes(x = PC, y = PER)) +
  geom_col(width = 0.5, color = "black") +
  xlab("Principal component") +
  ylab("Percentage of variation (%)") +
  theme_classic()

# Norm of the first eigenvector
norm(as.matrix(ev_1), "F")
## [1] 1
# Norm of the second eigenvector
norm(as.matrix(ev_2), "F")
## [1] 1

# Data frame with both eigenvectors
loads <- data.frame(
  VAR   = c("var_1", "var_2"),
  PC1 = ev_1, # First eigenvecor
  PC2 = ev_2  # Second eigenvectors
)
loads
##     VAR        PC1        PC2
## 1 var_1 -0.8134113  0.5816890
## 2 var_2 -0.5816890 -0.8134113
?solve
# Inverse of eigenvectors matrix
inv_evec <- solve(e_vec) 
# Change the basis of the original data 
data_set_3 <- data_set_2 %*% inv_evec
# Scatter showing the rotation 
ggplot(data.frame(data_set_3), aes(X1, X2)) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = 0, size = .5) +
  geom_hline(yintercept = 0, size = .5) +
  xlab("PC1 (78.8%)") +
  ylab("PC2 (21.2%)") +
  theme_classic()

library(ggpubr)
# Scatter plot with the centered data 
plot_data <- ggplot(data.frame(data_set_2), aes(x = varc_1, y = varc_2)) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = 0, size = .5) +
  geom_hline(yintercept = 0, size = .5) +
  ylim(c(-8, 8.5)) +
  ggtitle("Original Data") +
  theme_classic()
# Scatter plot with the rotated data
plot_rotation <- ggplot(data.frame(data_set_3), aes(X1, X2)) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = 0, size = .5) +
  geom_hline(yintercept = 0, size = .5) +
  xlab("PC1 (78.8%)") +
  ylab("PC2 (21.2%)") +
  ylim(c(-8, 8.5)) +
  ggtitle("Change of Basis to Eigenvectors") +
  theme_classic()
# Both graphs side by side
ggarrange(plot_data, plot_rotation)

# Data points just from PC 1
data_pc1 <- data.frame(v1 = data_set_3[,1], v2 = rep(0, nrow(data_set_3)))
# Scatter plot showing the projected points from PC1 (red points)
ggplot(data.frame(data_set_3), aes(X1, X2)) +
  geom_point(color = "blue", size = 2) +
  geom_point(data = data_pc1, aes(v1, v2), color = "red", size = 2) +
  geom_vline(xintercept = 0, size = .5) +
  geom_hline(yintercept = 0, size = .5) +
  xlab("PC1 (78.8%)") +
  ylab("PC2 (21.2%)") +
  ylim(c(-8, 8.5)) +
  theme_classic()
