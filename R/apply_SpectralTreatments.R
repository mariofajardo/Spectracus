#' Filtering methods for eSpectra objects
#' 
#' #' (NOT FULLY IMPLEMENTED) 
#'  
#' This method applies spectral pretreatments from a class \code{\link{eSpectra}} object 'A' to another \code{\link{eSpectra}} object 'B'
#' @rdname  apply_SpectralTreatments
#' @param eSpectra_A Object of class \code{\link{eSpectra}}.
#' @param eSpectra_B  
#' @author Mario Fajardo.
#' @examples 
#' \dontrun{
#' data("eSpectraExample")
#' filteredSpectra <- filter_eSpectra(eSpectraExample,c('SNV','C-hull'))
#' result <- apply_SpectralTreatments(eSpectraExample,filteredSpectra)
#' 
#' result
#' }
#' @exportMethod  apply_SpectralTreatments


setGeneric("apply_SpectralTreatments",
           function(eSpectra_A,eSpectra_B)
           {
             standardGeneric('apply_SpectralTreatments')
           }
)

setMethod(f = 'apply_SpectralTreatments',
          signature('eSpectra'),
          definition= function(eSpectra_A,eSpectra_B)
          {
            type <- gsub('(?<=\\s)[[:print:]]+','',eSpectra_B@Treatments, perl=T)
            TransformedSpectra <- filter_eSpectra(eSpectra_A,type = type)
            return(TransformedSpectra)
          }
) 



