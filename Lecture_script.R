library(codep)

# Load Doubs fish Data ----

data(Doubs)
species <- as.data.frame(Doubs.fish[-8,])
vars <- as.data.frame(cbind(Doubs.env[-8,],Doubs.geo[-8,]))

# Cottus gobio model
CHA.alt <- species$CHA + 1 #get rid of zeros
COGO_mod <- glm(CHA.alt ~ pH + flo + oxy, data = vars, family = Gamma(link = "log"))
summary(COGO_mod) 


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


