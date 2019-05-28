#' tidy method for eSpectra objects
#'
#'
#' @author Mario Fajardo
#' @param eSpectra A \code{\link{eSpectra}} object.
#' @return A ready to use longFormat data.frame suited for further manipulation by packages like dplyr, reshape2, ggplot2 and ggvis
#' 
#' @importFrom reshape2 melt
#'
#' @exportMethod tidy
#' @examples
#' \dontrun{
#' require(ggplot2)
#' data("eSpectraExample")
#'
#' TidySpec <- tidy(eSpectraExample)
#' 
#' head(TidySpec)
#'
#' ggplot(TidySpec)+geom_line(aes(x = nm, y = Reflectance, group = ID,colour=clay))
#'
#' }
#'
#' 
#' 
#' 
#' 

setGeneric('tidy',
           function(eSpectra)
             { 
             standardGeneric('tidy')
             }
           )

setMethod(f = 'tidy',
          signature="eSpectra",
          definition= function(eSpectra)
            {
              Specdata <- data.frame(ID=eSpectra@ID,eSpectra@Properties,eSpectra@Spectra)
              logic <- grepl('X',colnames(Specdata))
              SpecLong <- melt(Specdata,id.vars=colnames(Specdata)[!logic],variable.name = eSpectra@Units, value.name=eSpectra@Type)
              levels(SpecLong[,eSpectra@Units]) <- as.numeric(eSpectra@Bands)
              SpecLong[,eSpectra@Units] <- as.numeric(as.character(SpecLong[,eSpectra@Units]))
              SpecLong
            
          }
          )