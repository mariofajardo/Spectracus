#' points method for SoilSpectra objects
#' 
#' This method generates a plot.SoilSpectra for a Spectra Object
#' @rdname  points
#' @param SoilSpectra Object of class \code{\link{SoilSpectra}}.
#' @param ID numeric or character vector with ID of spectra to plot (default=1).
#' @author Mario Fajardo.
#' @exportMethod points
#' @examples  
#' \dontrun{
#' data("SoilSpectraExample")
#' plot(SoilSpectraExample)
#' points(SoilSpectraExample)
#' }


setMethod(f = 'points',
          signature(x='SoilSpectra'),
          definition= function(x,...)
            points(x@Bands,
                 x@Spectra,
                 ...)
          )