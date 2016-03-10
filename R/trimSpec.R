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

trimSpec <- function(spectra, wavlimits) {
  datawavs <- as.numeric(names(spectra))
  limits <- which(datawavs %in% wavlimits)
  kept_index <- seq(limits[1], limits[2], 1)
  trimmed_spectra <- spectra[, kept_index]
  kept_names <- datawavs[kept_index]
  colnames(trimmed_spectra) <- kept_names
  return(trimmed_spectra)}

