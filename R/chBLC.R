#' Continuum removal by removing the convex hull
#' 
#' This function applies a linearly interpolated convex hull to the spectra and returns the ratio of the deviation to the hull value
#' NOTE: This function will be deprecated in the next version, please use \code{\link{filter_spectra}} instead
#' @param spectra a matrix or data.frame with wavelengths as columns and spectra as rows
#' @export

chBLC <- function(spectra){
  warning('This function will be deprecated in the next version, please use filter_spectra() instead')
  interval <- seq_len(ncol(spectra))
  hull_spectra <- matrix(NA,ncol=ncol(spectra),nrow=nrow(spectra))
  for (i in seq_len(nrow(spectra))){
    tempSpect <- as.matrix(spectra[i,])
    data1 <- sortedXyData(interval, tempSpect)
    ## calculate convex hull
    c_hull <- chull(data1)
    ## get the appropriate region: the points of the polygon over the spectra

    # Create vector which wraps around
    c_hull <- c(c_hull, c_hull)
    # remove all points before the first one.
    c_hull <- c_hull[which.min(c_hull):length(c_hull)]
    # Go until the first end
    c_hull <- c_hull[1:which.max(c_hull)]
        
    ## calculate linear approximation between hull points
    linear_approx <- approx(data1[c_hull,], xout = interval, method = 'linear', ties = 'mean')
    ## calculate the deviation from the convex hull
    hull_spectra[i,] <- ( linear_approx[[2]] - tempSpect )/linear_approx[[2]]}
  colnames(hull_spectra) <- colnames(spectra)
  return(hull_spectra)}