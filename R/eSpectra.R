#' Formal S4 class eSpectra
#' 
#' @slot Meta A character vector with dataset details
#' @slot Instrument A character vector specifying de name of the Instrument where the spectram was obtained 
#' @slot Spectra A numeric matrix which contains absorbance/reflectance soil information 
#' @slot Bands A character vector specifying Band names for the spectra
#' @slot Units A character vector with the units of each Band.
#' @slot Id A character vector specyfing each spectrum Id.
#' @slot RowsAreSpectra  Logical 
#' @slot Type character vector specifying if the information is expressed in Absorbance, Reflectance, Energy or others (Normalized, EPO, etc.) .  
#' @slot Treatment Spectral Treatments that had being applied to the spectra
#' @name eSpectra
#' @rdname eSpectra
#' @exportClass eSpectra
#' @author Mario Fajardo, Brendan Malone, Budiman Minasny and Edward Jones.

setClass(Class = 'eSpectra',
         slots=c(Meta='character',
                 Instrument='character',
                 Spectra='matrix',
                 Bands='character',
                 Units='character',
                 RowsAreSpectra='logical',
                 Type='character',
                 ID='character',
                 Properties='data.frame',
                 Treatments=c('ANY'))
)



 