#' points method for eSpectra objects
#' 
#' This method generates a plot.eSpectra for a Spectra Object
#' @rdname  points
#' @param eSpectra Object of class \code{\link{eSpectra}}.
#' @param ID numeric or character vector with ID of spectra to plot (default=1).
#' @author Mario Fajardo.
#' @exportMethod points
#' @examples  
#' \dontrun{
#' data("eSpectraExample")
#' plot(eSpectraExample)
#' points(eSpectraExample)
#' }


setMethod(f = 'points',
          signature(x='eSpectra'),
          definition= function(x,...)
            points(x@Bands,
                 x@Spectra,
                 ...)
          )