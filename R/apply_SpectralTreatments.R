#' Filtering methods for SoilSpectra objects
#' 
#' This method applies spectral pretreatments from a class \code{\link{SoilSpectra}} object 'A' to another \code{\link{SoilSpectra}} object 'B'
#' @rdname  apply_SpectralTreatments
#' @param SoilSpectra_A Object of class \code{\link{SoilSpectra}}.
#' @param SoilSpectra_B  
#' @author Mario Fajardo.
#' @examples 
#' \dontrun{
#' data("SoilSpectraExample")
#' filteredSpectra <- filter_SoilSpectra(SoilSpectraExample,c('SNV','C-hull'))
#' result <- apply_SpectralTreatments(SoilSpectraExample,filteredSpectra)
#' 
#' result
#' }
#' @exportMethod  apply_SpectralTreatments


setGeneric("apply_SpectralTreatments",
           function(SoilSpectra_A,SoilSpectra_B)
           {
             standardGeneric('apply_SpectralTreatments')
           }
)

setMethod(f = 'apply_SpectralTreatments',
          signature('SoilSpectra'),
          definition= function(SoilSpectra_A,SoilSpectra_B)
          {
            type <- gsub('(?<=\\s)[[:print:]]+','',SoilSpectra_B@Treatments, perl=T)
            TransformedSpectra <- filter_SoilSpectra(SoilSpectra_A,type = type)
            return(TransformedSpectra)
          }
) 



