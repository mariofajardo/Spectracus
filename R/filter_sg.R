#' Sativsky-Golay filtering of spectra
#' 
#' This function performs a Sativsky Golay filter using \code{\link{sgolayfilt}} on matrices in which each row is a separate spectrum.
#' @param spectra matrix where each row is a spectrum and each column is a separate wavelength. It is important that these wavelengths be at least very nearly evenly spaced.
#' @param n,p,m parameters controlling \code{\link{sgolayfilt}} window size, polynomial order and derivative order respectively.
#' @note This function will be deprecated in the next version, please use \code{\link{filter_spectra}} instead
#' @importFrom plyr aaply
#' @importFrom signal sgolayfilt
#' @export filter_sg
#' @references Savitzky, Abraham and Golay, M. J. E. (1964) Smoothing and Differentiation of Data by Simplified Least Squares Procedures, Analytical Chemistry, 36, pp. 1627-1629.

filter_sg <-
function(spectra, n=11, p=2, m=1){
  warning('This function will be deprecated in the next version, please use filter_spectra() instead')
  spectra <- as.matrix(spectra)
  ## run filter
  sg <- aaply(spectra, 1, sgolayfilt, n = n, p = p, m = m)
  ##
  ## arrange appropriately if a single sample
  if(nrow(spectra) == 1){sg <- matrix(sg, dim(spectra))}
  ## return data frame
  colnames(sg) <- colnames(spectra)
  return(sg)
}
