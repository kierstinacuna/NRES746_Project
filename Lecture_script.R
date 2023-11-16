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


# Dissimilarity Measures ----

## Types of Distances Coefficients ----

### Euclidean Distance ----

# Simulate species abundance by site data set 
Y.hmm <- data.frame(y1 = c(0, 0, 1), y2 = c(4, 1, 0),
                    y3 = c(8, 1, 0), y4 = c(6, 0, 4))

# Calculate Euclidean distance using the dist() function
(Y.hmm.DistEu <- dist(x = Y.hmm, method = "euclidean"))

# Check the data type of the dist() function output
class(Y.hmm.DistEu)

# Coerce dist data type into a matrix
(eu_dist_matrix <- as.matrix(Y.hmm.DistEu))

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

# Unconstrained Ordination ----

## Principal Component Analysis ----

### Breaking Down a PCA ----

# Load the datasets package
library(datasets)

# Load the varechem dataset
data(varechem)

# Select data
(data <- varechem[, 1:2])

# Standardize the data to have mean zero and unit variance:
data_std <- scale(data)

# Compute the Covariance Matrix
cov_matrix <- cov(data_std)

# Perform the Eigendecomposition of the covariance matrix
eigen_decomp <- eigen(cov_matrix)
Eigenvalues <- eigen_decomp$values
Eigenvectors <- eigen_decomp$vectors

# Project the standardized data onto the Eigenspace
F_PrComps <- data_std %*% Eigenvectors
head(F_PrComps)

eig_vec = Eigenvectors
eig_val = Eigenvalues
std_data = data_std

# Plot the two principal components
plot_2pc <- function(eig_vec, eig_val, std_data) {
  eig_vec. <- as.data.frame(eig_vec)
  row.names(eig_vec.) <- c("P", "N")
  colnames(eig_vec.) <- c("PC1", "PC2")
  eig_vec.[, 1] <- eig_vec.[, 1]*-1
  Y_std <- as.data.frame(std_data)
  
  op <- par(mfrow = c(2, 1),     # 2x2 layout
            oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
            mar = c(1, 1, 0, 0), # space for one row of text at ticks and to separate plots
            mgp = c(2, 1, 0)    # axis label at 2 rows distance, tick labels at 1 row
  )       
  
  plot(N ~ P, 
       col = as.factor(rownames(Y_std)),
       pch = 19, 
       xlim=c(-2.2, 2.2), 
       ylim = c(-2.2,2.2), 
       data = as.data.frame(Y_std))
  
  abline(v=0 , h=0, 
         col = "dark gray")
  
  #Overlap pertinent evectors
  
  abline(0, 
         eig_vec[2, 1]/eig_vec[1, 1],
         col='purple')
  
  arrows(x0 = 0, 
         y0 = 0, 
         x1 = eig_val[1]*eig_vec[1, 1],
         y1 = eig_val[1]*eig_vec[2, 1],
         col = "purple",
         lwd = 2)
  
  # Plot the lines from first evector to points
  
  line1 <- c(0, 
             eig_vec[2, 1]/eig_vec[1, 1])
  
  perp.segment.coord <- function(x0, y0, line1){
    a <- line1[1]  #intercept
    b <- line1[2]  #slope
    x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
    y1 <- a + b * x1
    list(x0 = x0, y0 = y0, 
         x1 = x1, y1 = y1)
  }
  
  ss <- perp.segment.coord(Y_std$P, 
                           Y_std$N, 
                           line1)
  # do.call(segments, ss)
  # which is the same as:
  
  segments(x0 = ss$x0, 
           x1 = ss$x1, 
           y0 = ss$y0, 
           y1 = ss$y1, 
           col = 'purple')
  
  points(N ~ P, 
         col = as.factor(rownames(Y_std)), 
         pch = 19,
         data = Y_std)
  with(Y_std,
       text(N ~ P, 
            labels = as.factor(rownames(Y_std)),
            pos = 1, 
            cex=1.4))
  
  
  plot(N ~ P, 
       col = as.factor(rownames(Y_std)),
       pch = 19, 
       xlim=c(-2.2, 2.2), 
       ylim = c(-2.2,2.2), 
       data = as.data.frame(Y_std))
  
  abline(v=0 , h=0, 
         col = "dark gray")
  
  #Overlap pertinent evectors
  
  abline(0, 
         eig_vec[2, 1]/eig_vec[1, 1],
         col='purple')
  abline(0, 
         eig_vec[1, 2]/eig_vec[2, 2],
         col='orange')
  
  arrows(x0 = 0, 
        y0 = 0, 
        x1 = eig_val[1]*eig_vec[1, 1],
         y1 = eig_val[1]*eig_vec[2, 1],
         col = "purple",
         lwd = 2)
  
  arrows(x0 = 0, 
         y0 = 0, 
         x1 = eig_val[2]*eig_vec[1,2], 
         y1 = eig_val[2]*eig_vec[2, 2],
         col = "orange", 
         lwd = 2)
  
  
  line2 <- c(0, 
             eig_vec[1, 2]/eig_vec[1, 1])
  
  perp.segment.coord <- function(x0, y0, line2){
    a <- line2[1]  #intercept
    b <- line2[2]  #slope
    x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
    y1 <- a + b * x1
    list(x0 = x0, y0 = y0, 
         x1 = x1, y1 = y1)
  }
  
  ss <- perp.segment.coord(Y_std$P, 
                           Y_std$N, 
                           line2)
  
  segments(x0 = ss$x0, 
           x1 = ss$x1, 
           y0 = ss$y0, 
           y1 = ss$y1, 
           col = 'orange')
  
  points(N ~ P, 
         col = as.factor(rownames(Y_std)), 
         pch = 19,
         data = Y_std)
  
  with(Y_std,
       text(N ~ P, 
            labels = as.factor(rownames(Y_std)),
            pos = 1, 
            cex=1.4)
  )
  
  title(xlab = "N",
        ylab = "P",
        outer = TRUE, line = 3)
}

plot_2pc(Eigenvectors, Eigenvalues, data_std)
