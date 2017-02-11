# Function_to_calculate_eigenvectors
# Alexey Larionov, 07Feb2017

# Notes:
#
# This script implements procedure described by Price et al 2006 (PMID: 16862161).
# Note that Patterson et al 2006 (PMID: 17194218) recommend slightly different normalisation. 
#
# The function takes a numeric matrix with genotypes (e.g. coded as 0, 1, 2, NA).
# The matrix should have cases in columns and variants in rows.
# The function outputs and "eigen" object, which is a list containing 
# - matrix of eigenvectors and 
# - vector of eigenvalues (see details in ?eigen)

normalise_and_calculate_eigenvectors <- function(x) {
  
  # --- Center and normalise variants (rows) --- #
  
  # Center by mean
  avg.rows <- apply(x, 1, mean, na.rm=TRUE)
  x.c <- x - avg.rows
  
  # Normalise by sqrt(p(1-p)) where p~"posterior estimate of unobserved allele frequency"
  # This is motivated by the fact that genetic drift per generation is proportional to this normalisation value (Patterson 2006)
  # Also this makes each column to have same variance
  # 
  p.fnc <- function(x) (1 + sum(x, na.rm=TRUE)) / (2 + 2 * sum(!is.na(x)))
  p <- apply(x, 1, p.fnc)
  eaf <- sqrt(p*(1-p))
  x.cn <- x.c/eaf
  
  # Substitute NAs to zeros
  0 -> x.cn[is.na(x)]
  
  # --- Calculate eigenvectors of covariance matrix of cases --- #
  
  cov.mx <- cov(x.cn)
  eig <- eigen(cov.mx) # eigenvectors in columns
  
  return(eig)
  
}

