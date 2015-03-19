#' Wavelet smoothing
#' 
#' @param spectra a matrix or data.frame where the columns are wavelengths and the rows, individual spectra
#' @param res level to be extracted from wavelet decomposition model
#' 
#' @author Brendan Malone
#' 
#' @importFrom wavethresh wd accessC.wd
#' @export

wavelet_smooth <-function(spectra, res){
  wave_spectra<- matrix(NA,ncol=2^res,nrow=nrow(spectra))
  for (i in 1:nrow(spectra)){
    wds<-wd(as.matrix(spectra[i,]),bc="symmetric",filter.number = 10, family = 'DaubExPhase', min.scale = 2)
    wave_spectra[i,]<- accessC.wd(wds, level=res)}
  colnames(wave_spectra) <- seq((404 + 0.5*(2048/(2^res))),2451, by=2048/(2^res))
  return(wave_spectra)
}