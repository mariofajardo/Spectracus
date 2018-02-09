#' print method for SoilSpectra objects
#' 
#' This method generates a description for a Spectra Object
#' @rdname  print
#' @param SoilSpectra Object of class \code{\link{SoilSpectra}}.
#' @author Mario Fajardo.
#' @exportMethod   print
#' @examples  
#' \dontrun{
#' data("SoilSpectraExample")
#' print(SoilSpectraExample)
#' }



setMethod(f = 'print',
          signature(x='SoilSpectra'),
          definition= function(x)
          {
            
            cat('\n',length(x@ID),' soil observations','\n\n',
                'Spectra values ranging from ~', round(as.numeric(x@Bands)[1]),' to ',round(tail(as.numeric(x@Bands))[6]),x@Units,'\n\n',
                'Data type : ',paste(x@Treatments,collapse = '+'),'\n\n',
                ncol(x@Properties),'associated attributes')
          }
) 
