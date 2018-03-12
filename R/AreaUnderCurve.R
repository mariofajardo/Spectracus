#' Spectra normalisation via continuum estimation
#' 
#' Generic Function for:
#' \itemize{
#' \item 1. Fitting continuum from spectral hull features
#' \item 2. Outputting continuum removed or normalised spectra
#' \item 3. Estimating normalised hull feature polygon vertices
#' }
#' 
#' This function is useful in studies for example where spectral features of phenomena need to be investigated
#' Selecting the convex hull features in the right order has been an ongoing teething issue.
#' 
#' @aliases AreaUnderCurve
#' @importFrom plyr adply
#' @importFrom plyr splat
#' @importFrom plyr llply
#' @importFrom munsell rgb2mnsl
#' @importFrom splancs areapl
#' 
#' 
#' 
#' @param Spectra A \code{\link{SoilSpectra}} object.
#' @param AUC A logical that is expressed as the Area Under the Curve (AUC) or Over the Curve is desired, default is set to TRUE i.e., AUC
#' @return Returns a list with 8 elements:
#' \itemize{
#' \item{hull}{ A \code{\link{SoilSpectra}} object with hull spectra (equivalent to \code{\link{filterSpectra}} output)}
#' \item{raw}{ A \code{\link{SoilSpectra}} object with the input spectra}
#' \item{polygons}{ Vertices for constructing the polygons of each normalised spectral feature}
#' \item{contSpectra} {A \code{\link{SoilSpectra}} object with the Continuum removed spectra}
#' \item{areas}{ areas of the polygons}
#' \item{depths}{ difference between the highest reflectance/absorbance and the lower of each polygon}
#' \item{slopes}{ slopes between the highest reflectance/absorbance and the lower of each polygon}
#' 
#' }
#' @examples 
#' \dontrun{
#' data("SoilSpectraExample")
#' #Select a specific area
#' SelectedArea <- SoilSpectraExample[122:123,1950:2220]
#' 
#' 
#' #Run the continuum fitting function for Area Under the Curve
#' spec.CR<- AreaUnderCurve(SelectedArea,AUC=T)
#' 
#' 
#' plot(SoilSpectraExample[122:123,])
#' points(spec.CR$raw[1],col='blue',type='l',lwd=2)
#' points(spec.CR$raw[2],col='red',type='l',lwd=2)
#' 
#' 
#' points(spec.CR$continuum[1],col='blue',type='l')
#' points(spec.CR$continuum[2],col='red',type='l')
#' 
#' 
#' 
#' #check Differences
#' plot(spec.CR$hull)
#' polygon(spec.CR$polygons$`122`,col='blue')
#' polygon(spec.CR$polygons$`123`,col='red')
#' 
#' spec.CR$areas
#' spec.CR$slopes
#' spec.CR$depths
#' 
#' 
#' 
#' 
#' 
#' #Run the continuum fitting function for Area Over the curve
#' SelectedArea <- SoilSpectraExample[122:123,2180:2230]
#' 
#' spec.CR<- AreaUnderCurve(SelectedArea,AUC=F)
#' 
#' 
#' plot(SoilSpectraExample[122:123,])
#' points(spec.CR$raw[1],col='blue',type='l',lwd=2)
#' points(spec.CR$raw[2],col='red',type='l',lwd=2)
#' 
#' 
#' points(spec.CR$continuum[1],col='blue',type='l')
#' points(spec.CR$continuum[2],col='red',type='l')
#' 
#' 
#' 
#' #check Differences
#' plot(spec.CR$hull)
#' polygon(spec.CR$polygons$`122`,col='blue')
#' polygon(spec.CR$polygons$`123`,col='red')
#' 
#' spec.CR$areas
#' spec.CR$slopes
#' spec.CR$depths
#' }
#' 
#' @exportMethod  AreaUnderCurve
#' @author Mario Fajardo and Brendan Malone



setGeneric("AreaUnderCurve",
           function(SoilSpectra,AUC=TRUE,...)
           {
             standardGeneric('AreaUnderCurve')
           }
)

setMethod(f = 'AreaUnderCurve',
          signature = 'SoilSpectra',
          definition= function(SoilSpectra,AUC=TRUE,...)
            {
            
            spectra <- SoilSpectra@Spectra
            interval <- SoilSpectra@Bands
            
            
            result <- apply(spectra,1,function(tempSpectrum) {
              
              WarnFlag <- NA
              if(any(tempSpectrum<0)){
              correction<-tempSpectrum[tempSpectrum<0][which.min(tempSpectrum[tempSpectrum<0]-0)]
              under_0 <- interval[tempSpectrum<0]
              tempSpectrum <- tempSpectrum+abs(correction)
              WarnFlag <- paste0('Please check your spectra, it had negative values')
              } 
              
              hull_spectra <- matrix(NA, ncol = length(tempSpectrum), nrow = 1)
              tempSpect <-  round(matrix(tempSpectrum,ncol = length(tempSpectrum)),digits = 4)
              
              tmpSortedData <- sortedXyData(interval,tempSpectrum)
              tmpSortedData$y <- round(tmpSortedData$y,digits = 4)
              # plot(tmpSortedData$x,tmpSortedData$y,type='b')
              c_hull <- chull(tmpSortedData)
              # points(tmpSortedData$x[c_hull],tmpSortedData$y[c_hull],col='red',cex=2)
              first_last <- data.frame(x=tmpSortedData[c(1,nrow(tmpSortedData)),1],y=tmpSortedData[c(1,nrow(tmpSortedData)),2])  
              linear_model <- lm(formula = y~x,data = first_last)
              
              if(!AUC){ #if AUC select those points higher than the linear approx between extremes
                is_up <- as.logical((tmpSortedData[c_hull,2]>=round(predict(linear_model,tmpSortedData[c_hull,]),digits = 4)))
                c_hull <- c_hull[is_up]
                }else{
                  is_down <- as.logical((tmpSortedData[c_hull,]<=round(predict(linear_model,tmpSortedData[c_hull,]),digits = 4))[,2])
                  c_hull <- c_hull[is_down]
                  }
              # points(tmpSortedData[c_hull,],col='blue',pch=16)
              # lines(tmpSortedData[c_hull,][order(tmpSortedData[c_hull,1]),],col='blue',pch=16)
              
              linear_approx <- approx(tmpSortedData[c_hull, ], xout = interval,
                                      method = "linear", ties = "mean")
              
              linear_approx$y <- round(linear_approx$y,digits = 4)
              if (!AUC) {
                hull_spectra[1, ] <- 1-((linear_approx$y - tempSpectrum)/linear_approx$y)
                hull_spectra[, c(1, ncol(hull_spectra))] <- 1.0001 #in order to make the extreme points unique
                hull_data <- data.frame(wave=interval,value=t(hull_spectra))
                closing_line <- data.frame(wave=hull_data$wave, value=rep(1.0001, each = length(hull_data$wave)))
                
                }else{
                  
                hull_spectra[1, ] <- (tempSpect-linear_approx[[2]])/linear_approx[[2]]
                hull_spectra[, c(1, ncol(hull_spectra))] <- -0.0001 #in order to make the extreme points unique
                hull_data <- data.frame(wave=interval,value=t(hull_spectra))
                hull_data$value[is.nan(hull_data$value)] <- 0
                closing_line <- data.frame(wave=hull_data$wave, value=rep(-0.0001, each = length(hull_data$wave)))
                }
              # plot(hull_data,type='b')
              # points(closing_line,col='green',type='b')
              
              hull_polygon <- as.matrix(rbind(hull_data, closing_line))
              # polygon(hull_polygon,col='red')
              
              if(length(c_hull)<2){
                area_feature <- NA 
                slope_feature <- linear_model$coefficients[2]
                }else{
                  area_feature <- areapl(hull_polygon)
                  slope_feature <- linear_model$coefficients[2]
                }
              
              if(!AUC){
                min_raw <- tempSpect[which(hull_spectra==min(hull_spectra))[1]]
                min_continuum <- hull_data$value[which(hull_spectra==min(hull_spectra))[1]]
                if(min_raw==0&min_raw==0){ #check
                  depth_feature <- 1
                  }else{
                    depth_feature <- 1-(min_raw/min_continuum)
                  }
                
                }else{
                  min_raw <- tempSpect[which(hull_spectra==min(hull_spectra))[1]]
                  min_continuum <- hull_data$value[which(hull_spectra==min(hull_spectra))[1]]
                  depth_feature <- tempSpect[which(hull_spectra==min(hull_spectra))[1]]/hull_data$value[which(hull_spectra==min(hull_spectra))[1]]
                }
              
              
              retval <- list(wave = hull_data$wave,
                             c.hull = hull_spectra,
                             raw.spec = tmpSortedData$y,
                             continuum = linear_approx[[2]],
                             polygon = hull_polygon,
                             area=area_feature,
                             depth=depth_feature,
                             slope=slope_feature,
                             warnFlag=WarnFlag)
              })
            names(result) <- SoilSpectra@ID

            hullSpectra <- initialize(SoilSpectra,
                                      Meta=x@Meta,
                                      Instrument=SoilSpectra@Instrument,
                                      Spectra=do.call(rbind,lapply(result,function(x) x$c.hull)),
                                      Bands=SoilSpectra@Bands,
                                      Units=SoilSpectra@Units,
                                      RowsAreSpectra=SoilSpectra@RowsAreSpectra,
                                      Type=SoilSpectra@Type,
                                      ID=SoilSpectra@ID,
                                      Properties=SoilSpectra@Properties,
                                      Treatments=SoilSpectra@Treatments)
            
            rawSpectra <- initialize(SoilSpectra,
                                     Meta=x@Meta,
                                      Instrument=SoilSpectra@Instrument,
                                      Spectra=do.call(rbind,lapply(result,function(x) x$raw.spec)),
                                      Bands=SoilSpectra@Bands,
                                      Units=SoilSpectra@Units,
                                      RowsAreSpectra=SoilSpectra@RowsAreSpectra,
                                      Type=SoilSpectra@Type,
                                      ID=SoilSpectra@ID,
                                      Properties=SoilSpectra@Properties,
                                      Treatments=SoilSpectra@Treatments)
            
            contSpectra <- initialize(SoilSpectra,
                                      Meta=x@Meta,
                                     Instrument=SoilSpectra@Instrument,
                                     Spectra=do.call(rbind,lapply(result,function(x) x$continuum)),
                                     Bands=SoilSpectra@Bands,
                                     Units=SoilSpectra@Units,
                                     RowsAreSpectra=SoilSpectra@RowsAreSpectra,
                                     Type=SoilSpectra@Type,
                                     ID=SoilSpectra@ID,
                                     Properties=SoilSpectra@Properties,
                                     Treatments=SoilSpectra@Treatments)
            
            polygons <- lapply(result,function(x) x$polygon)
            
            areas <- do.call(c,lapply(result,function(x) x$area))
            
            slopes <- do.call(c,lapply(result,function(x) x$slope))
            names(slopes) <- SoilSpectra@ID
            
            depths <- do.call(c,lapply(result,function(x) x$depth))
            
            FinalList <-list(hull=hullSpectra,
                             raw=rawSpectra,
                             continuum=contSpectra,
                             polygons=polygons,
                             areas=areas,
                             slopes=slopes,
                             depths=depths) 
            
            return(FinalList)
          }
)

