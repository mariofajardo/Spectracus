#' Wavelet smoothing
#' 
#' @param spectra a matrix or data.frame where the columns are wavelengths and the rows, individual spectra
#' @param res level (an integer >1) to be extracted from wavelet decomposition model
#' 
#' @author Brendan Malone
#' 
#' @importFrom wavethresh wd accessC.wd
#' @export

# FUNCTION: wavelet smoothing 
wavelet_smooth <- function(spectra, res) { 
  nm2<- 2^c(1:100)
  vs<- ncol(spectra)
  if (sum(nm2 == vs) != 1) {stop("Error: Number of columns in spectra table needs to equal 2^x")} else {wave_spectra <- matrix(NA, ncol = 2^res, nrow = nrow(spectra)) 
  for (i in 1:nrow(spectra)){ 
    wds <- wd(as.matrix(spectra[i, ]), bc = "symmetric", filter.number = 10, family = "DaubExPhase", min.scale = 2) 
    wave_spectra[i, ] <- accessC.wd(wds, level = res)} 
  wave_spectra<- as.data.frame(wave_spectra)
  colnames(wave_spectra) <- seq((as.numeric(names(spectra)[1]) + 0.5 * (ncol(spectra)/(2^res))), as.numeric(names(spectra)[length(spectra)]), by = ncol(spectra)/(2^res)) }
  return(wave_spectra)}