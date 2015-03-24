#' Trim regular spectra
#' NOTE: Need to make this more general
#' Convenience wrapper for \link{strip_spectra}, for wavelengths 350:2500
#' 
#' @param spectra a matrix of spectra where each row represents a spectrum and each column a wavelength
#' @param wavlimits minimum and maximum wavelengths desired. Must be present in \code{datawavs}
#' @param datawavs a numeric or integer vector of all wavelengths
#' @param ... Further arguments to be pass to \link{strip_spectra}
#' 
#' @export

trimSpec <- function(spectra, wavlimits, datawavs=350:2500, ...){

  strip_spectra(spectra=spectra, wavlimits=wavlimits, datawavs=datawavs, ...)
  
}

# trimSpec<- function (spectra, wavlimits, wavelength=350:2500){
#   datawavs <- as.numeric(wavelength)  # make a long sequence of all wavelength labels in the spectrum
#   limits <- which(datawavs %in% wavlimits)  #Which column labels match up with the wavelength bounds
#   kept_index <- seq(limits[1], limits[2], 1)
#   trimmed_spectra <- spectra[,kept_index] # The trimmed spectra
#   kept_names <- datawavs[kept_index]
#   colnames(trimmed_spectra) <- kept_names
#   return(trimmed_spectra)}