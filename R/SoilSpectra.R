#' Formal S4 class SoilSpectra
#' 
#' @slot Instrument A character vector specifying de name of the Instrument where the spectram was obtained 
#' @slot Spectra A numeric matrix which contains absorbance/reflectance soil information 
#' @slot Wavelength A character vector specifying wavlength names
#' @slot Range A character vector section of the electromagnetic spectrum used e.g., Vis-NIR, NIR, MIR, etc.
#' @slot Wavenumber A character vector specifying Wavenumber names
#' @slot Id A character vector specyfing each spectrum Id.
#' @slot RowsAreSpectra  Logical 
#' @slot Type character vector specifying if the information is expressed in Absorbance, Reflectance or Modified (Normalized, EPO, etc.) .  
#' @slot Treatment Spectral Treatments that had being applied to the spectra
#' @name SoilSpectra
#' @rdname SoilSpectra
#' @exportClass SoilSpectra
#' @author Mario Fajardo, Brendan Malone, Budiman Minasny and Edward Jones.

setClass(Class = 'SoilSpectra',
         slots=c(Instrument='character',
                 Spectra=c('array'),
                 Wavelength='numeric',
                 Range='character',
                 Wavenumber='numeric',
                 RowsAreSpectra='logical',
                 Type='character',
                 ID='character',
                 Properties='data.frame',
                 Treatments=c('ANY'))
)



 