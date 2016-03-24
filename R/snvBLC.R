#' Standard Normal Variate (SNV) transformation
#' @param spectra a matrix or data.frame with wavelengths as columns and spectra as rows
#' @note This function will be deprecated in the next version, please use \code{\link{filter_spectra}} instead
#' @export

snvBLC<- function(spectra){
  warning('This function will be deprecated in the next version, please use filter_spectra() instead')
  spectra<-as.matrix(spectra)
  snvMat<-(spectra - rowMeans(spectra))/apply(spectra,1,sd)
  attributes(snvMat) <- attributes(spectra)
  snvMat
}

