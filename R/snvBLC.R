#' Standard Normal Variate (SNV) transformation
#' 
#' @param spectra a matrix or data.frame with wavelengths as columns and spectra as rows
#' @export

snvBLC<- function(spectra){
  spectra<-as.matrix(spectra)
  snvMat<-(spectra - rowMeans(spectra))/apply(spectra,1,sd)
  attributes(snvMat) <- attributes(spectra)
  snvMat
}