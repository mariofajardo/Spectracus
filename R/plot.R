#' plot method for SoilSpectra objects
#' 
#' This method generates a plot.SoilSpectra for a Spectra Object
#' @rdname  plot
#' @param SoilSpectra Object of class \code{\link{SoilSpectra}}.
#' @param ID numeric or character vector with ID of spectra to plot (default=1).
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
          {if(nrow(x@Spectra)==1){
            plot(x@Bands,
                x@Spectra,
                ylab=x@Type,
                xlab=paste0('Units: ',x@Units),
                type='l',
                main=paste0('Sample ',x@ID),
                ...)
          }else{
            plotMatrix <- t(x@Spectra)
            
            matplot(x@Bands,
                    plotMatrix,
                    ylab=x@Type,
                    xlab=paste0('Units: ',x@Units),
                    type='l',
                    ylim=range(x@Spectra),
                    ...)
            
          }
           
          }
) 

