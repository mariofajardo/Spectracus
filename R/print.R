#' print method for eSpectra objects
#' 
#' This method generates a description for a Spectra Object
#' @rdname  print
#' @param eSpectra Object of class \code{\link{eSpectra}}.
#' @author Mario Fajardo.
#' @exportMethod   print
#' @examples  
#' \dontrun{
#' data("eSpectraExample")
#' print(eSpectraExample)
#' }



setMethod(f = 'print',
          signature(x='eSpectra'),
          definition= function(x)
          {
            
            cat('\n',x@Meta,'\n\n',length(x@ID),' soil observations','\n\n',
                'Spectra values ranging from ~', round(as.numeric(x@Bands)[1]),' to ',round(tail(as.numeric(x@Bands))[6]),x@Units,'\n\n',
                'Data Treatments : ',paste(x@Treatments,collapse = '+'),'\n\n',
                ncol(x@Properties),'associated attributes')
          }
) 
