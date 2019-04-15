#' Perform Empirical Parameter Orthogonalisation
#' 
#' (NOT YET IMPLEMENTED) 
#'  
#' This function takes in a difference matrix (the difference between wet 
#' and dry soil spectra) and returns an EPO-transformed projection matrix.
#' 
#' @author Budiman Minasny
#' 
#' @param D matrix. Difference matrix between wet and dry soil spectra.
#' @param npc number of components to use from the Singular Value Decomposition.
#' 
#' @return A square matrix with dimensions equal to the number of columns of the difference matrix.
#' 
#' @references Roger, J.M., Chauchard, F., Bellon-Maurel, V., 2003. EPO-PLS external parameter orthogonalisation of PLS application to temperature-independent measurement of sugar content of intact fruits. Chemometrics and Intelligent Laboratory Systems 66, 191-204.
#'
#' Minasny, B., McBratney, A. B., Bellon-Maurel, V., Roger,J.-M., Gobrecht, A., Ferrand, L., Joalland, S., 2011. Removing the effect of soil moisture from NIR diffuse reflectance spectra for the prediction of soil organic carbon. Geoderma 167-168, 118-124.
#' 
#' @export

epo <- function(D, npc){
  # D is the difference matrix (spectra difference between wet and dry)
  # npc is the number of components to use
  # return: P: the projection matrix
  D <- as.matrix(D)
  n <- nrow(D)
  p <- ncol(D)
  
  dtd <- t(D) %*% (D)
  s <- svd(dtd) # singular value decomposition of the D (n x n) matrix
  ld <- s$v[, 1:npc]         # extract the no. factors
  P <- diag(p) - ld %*% t(ld)  # projection matrix
  
  P
}
