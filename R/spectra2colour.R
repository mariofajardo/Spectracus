#' Method that extracts RGB and Munsell colour from \code{\link{SoilSpectra}} objects.
#' 
#' Converts spectra reflectance into RGB and Munsell colours
#' @importFrom plyr adply
#' @importFrom plyr splat
#' @importFrom plyr llply
#' @importFrom munsell rgb2mnsl
#' 
#' @param SoilSpectra A \code{\link{SoilSpectra}} object.
#' 
#' @exportMethod spectra2colour
#' @author Michael Nelson and Mario Fajardo
#' @examples 
#' \dontrun{
#' data("SoilSpectraExample")
#' 
#' spectra2colour(SoilSpectraExample)
#' 
#' spectra2colour(SoilSpectraExample[23]) 
#' }
setGeneric("spectra2colour",
           function(SoilSpectra,...)
           {
             standardGeneric('spectra2colour')
           }
)

setMethod(f = 'spectra2colour',
          signature = 'SoilSpectra',
          definition= function(SoilSpectra,...){
            
            #Internal Functions
            spectra_to_RGB <- function(.spectra, 
                                       all_wavelengths,
                                       rgb_list = list(blue = list(.interval=450:520), 
                                                       red = list(.interval=600:690), 
                                                       green=list(.interval=520:600)), 
                                       ...) {
              ## a function to return the average values in the
              ## red, green and blue sections of the spectrum
              ## would work on any intervals
              ##
              ## get the appropriate indices
              interval_list <- llply(rgb_list, splat(in_interval), .all=all_wavelengths)
              ##
              ## get the average in these subsets
              rgb_values <- lapply(interval_list, mean_interval, .data=.spectra)
              ##
              ## convert to colour
              colour <- with(rgb_values, rgb(red, green, blue))
              ##
              ## return data frame
              with(rgb_values, data.frame(red, green, blue, colour))
            }
            
            in_interval <- function(.all, .interval,...){
              ## index for an interval
              ## a function to subset a particular waveband interval
              .in_interval = .all %in% .interval
            }
            
            mean_interval <- function(.data, .index){
              ## returns the mean for given indices
              mean(.data[.index])
            }
            
            ##End of internal functions
            
            spectra <- SoilSpectra@Spectra
  
  
  
  
  ## find r,g,b colour
  rgb_colours <- adply(spectra, 1, spectra_to_RGB, all_wavelengths = SoilSpectra@Bands,.id = NULL)
  ##
  ## get munsell colour
  munsell_colours <- splat(function(red,green,blue, ...){rgb2mnsl(R=red,G=green,B=blue)})(rgb_colours)
  ##
  ## return
  soil_colours <- data.frame(ID=SoilSpectra@ID,rgb_colours, munsell = munsell_colours)
  
  return(soil_colours)
          }
)




