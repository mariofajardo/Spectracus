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
#' \item{ContRemoved}{The normalised or continuum removed spectra}
#' \item{InputSpec}{The input spectra}
#' \item{continuum}{The fitted continuum}
#' \item{polygon}{Vertices for constructing the polygon of the normalised spectral feature}
#' \item{area}{area of the polygon}
#' \item{depth}{difference between the highest reflectance/absorbance and the lower}
#' \item{slope}{slope between the highest reflectance/absorbance and the lower}
#' 
#' }
#' @examples 
#' \dontrun{
#' data("SoilSpectraExample")
#' #Select a specific area
#' SelectedArea <- SoilSpectraExample[123,2180:2230]
#' 
#' plot(SelectedArea)
#' 
#' #Run the continuum fitting function
#' spec.CR<- AreaUnderCurve(SelectedArea,AUC=F)
#' 
#' 
#' plot(filter_SoilSpectra(SoilSpectraExample[123,2000:2500],'S-Golay'))
#' points(SelectedArea@Wavelength,spec.CR$`123`$raw.spec,type="l", xlab="wavelength", ylab="reflectance",col='red',lwd=5)
#' lines(SelectedArea@Wavelength,spec.CR$`123`$continuum,col="red",lwd=5)
#' 
#' #Plot continuum removed spectrum
#' plot(SelectedArea@Wavelength,spec.CR$`123`$c.hull,type="l", xlab="wavelength", ylab="reflectance")
#' polygon(spec.CR$`123`$polygon,col='red')
#' 
#' 
#' 
#' SelectedArea <- SoilSpectraExample[123,2000:2180]
#' 
#' plot(SelectedArea)
#' 
#' #Run the continuum fitting function
#' spec.CR<- AreaUnderCurve(SelectedArea,AUC=T)
#' 
#' 
#' plot(filter_SoilSpectra(SoilSpectraExample[123,2000:2500],'S-Golay'))
#' points(SelectedArea@Wavelength,spec.CR$`123`$raw.spec,type="l", xlab="wavelength", ylab="reflectance",col='red',lwd=5)
#' lines(SelectedArea@Wavelength,spec.CR$`123`$continuum,col="red",lwd=5)
#' 
#' #Plot continuum removed spectrum
#' plot(SelectedArea@Wavelength,spec.CR$`123`$c.hull,type="l", xlab="wavelength", ylab="reflectance")
#' polygon(spec.CR$`123`$polygon,col='red')
#' }
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
          definition= function(SoilSpectra='SoilSpectra',AUC=TRUE,...)
            {
            
            spectra <- SoilSpectra@Spectra
            interval <- SoilSpectra@Wavelength
            
            
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
            return(result)
          }
)
