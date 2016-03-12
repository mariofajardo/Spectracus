#' Reduce column dimensions
#' 
#' This function reduces the column dimensions by averaging values inside a window of a specified size
#' 
#' 
#' @param spectra matrix where each row is a spectrum and each column a wavelength
#' @param window window size over which the spectra will be averaged
#' @examples
#'#load data
#'data("soilSpecDat")
#'
#'#subset spectra
#'spectra <- soilSpecDat[,grepl('[X]',colnames(soilSpecDat))]
#'colnames(spectra) <- substr(colnames(spectra),2,5)
#'wavelengths <- as.numeric(colnames(spectra))
#'
#'#plot spectrum
#'plot(wavelengths,spectra[1,],type='b')
#'
#'spectra_2 <- compSpec(spectra,9)
#'wavelength_2 <- colnames(spectra_2)
#'
#'#plot processed spectrum
#'points(wavelength_2,spectra_2[1,],cex=2,col='red',pch=21,bg='blue') 
#'
#'
#' @author Brendan Malone
#' @export

compSpec <- function(spectra, window) {
  if(ncol(spectra)%%window != 0) 
  {stop("Error: Pick a more compatable window size!")} else {compMat <- matrix(NA, ncol = (ncol(spectra))/window, nrow = nrow(spectra))
  cc <- 1
  for (i in 1:ncol(compMat)) {
    compMat[, i] <- rowMeans(spectra[, cc:(cc + (window - 1))])
    cc <- cc + window}
  colab = seq(as.numeric(names(spectra)[1]), as.numeric(names(spectra)[length(spectra)]), by = window)
  compMat<- as.data.frame(compMat)
  colnames(compMat) <- colab}
  return(compMat)
}