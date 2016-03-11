#' Reduce column dimensions
#' 
#' This function reduces the column dimensions by averaging values inside a window of a specified size
#' 
#' 
#' @param spectra matrix where each row is a spectrum and each column a wavelength
#' @param window window size over which the spectra will be averaged
#' @examples
#' observed <- c(1,5,6,8,9,20)
#' predicted <- c(2,5,6,7,9,22)
#' goof(observed,predicted)
#' goof(observed,predicted,coefficient = c('R2','RMSE','RPD'),cex=3,pch=19,col='blue')
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