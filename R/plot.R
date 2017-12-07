#' plot method for SoilSpectra objects
#' 
#' This method generates a plot for a Spectra Object
#' @rdname  plot
#' @param SoilSpectra Object of class \code{\link{SoilSpectra}}.
#' @author Mario Fajardo.
#' @exportMethod plot
#' @examples  
#' \dontrun{
#' data("SoilSpectraExample")
#' plot(SoilSpectraExample)
#' plot(SoilSpectraExample,type='l')
#' }


setMethod(f = 'plot',
          signature(x='SoilSpectra'),
          definition= function(x,...)
          {
            plot(x@Wavelength,
                x@Spectra[1,],
                ylab=x@Type,
                xlab=slotNames(x)[3],
                ...)
          }
) 

