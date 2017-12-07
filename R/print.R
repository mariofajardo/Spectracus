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
                x@Type,'values ranging from', x@Wavelength[1],' to ',tail(x@Wavelength)[6],'nm in the ',x@Range,'region of the electromagnetic spectrum','\n\n',
                'Data type : ',x@Treatments)
          }
) 
