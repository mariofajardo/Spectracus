#' show method for SoilSpectra objects
#' 
#' This method generates a description for a Spectra Object
#' @rdname  show
#' @param SoilSpectra Object of class \code{\link{SoilSpectra}}.
#' @author Mario Fajardo.
#' @exportMethod  show


setMethod(f = 'show',
          signature('SoilSpectra'),
          definition= function(object)
          {
            print(object)
          }
) 

