#' show method for eSpectra objects
#' 
#' This method generates a description for a eSpectra Object
#' @rdname  show
#' @param eSpectra Object of class \code{\link{eSpectra}}.
#' @author Mario Fajardo.
#' @exportMethod  show


setMethod(f = 'show',
          signature('eSpectra'),
          definition= function(object)
          {
            print(object)
          }
) 

