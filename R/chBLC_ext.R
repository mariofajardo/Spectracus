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
#' @aliases chBLC_ext
#' @importFrom plyr adply
#' @importFrom plyr splat
#' @importFrom plyr llply
#' @importFrom munsell rgb2mnsl
#' 
#' 
#' 
#' @param spectra A single spectrum of class \code{data.frame}. Column values represent spectral values per wavelength unit.
#' @param lower numeric; The lower wavelength bound of the spectral feature or region of interest.
#' @param upper numeric; The upper wavelength bound of the spectral feature or region of interest.
#' @param AUC A logical that is expressed as the Area Under the Curve (AUC) or Over the Curve is desired, default is set to TRUE i.e., AUC}
#' 
#' @return Returns a list with 8 elements:
#' \itemize{
#' \item{wave}{ The sequence of spectral wavelengths of the region of interest.}
#' \item{c.hull}{ The normalised or continuum removed spectrum}
#' \item{raw.spec}{ The raw spectrum}
#' \item{continuum}{ The fitted continuum}
#' \item{polygon}{ Vertices for constructing the polygon of the normalised spectral feature}
#' \item{area}{ NEEDS A GOOD DESCRIPTION AND REFERENCES HOPEFULLY}
#' \item{depth}{ NEEDS A GOOD DESCRIPTION AND REFERENCES HOPEFULLY}
#' \item{slope}{ NEEDS A GOOD DESCRIPTION AND REFERENCES HOPEFULLY}
#' 
#' }
#' @note {This function works similarly to the other functions in this package where a continuum is fitted to data from which raw spectra values are removed from it. However, this function provides extra useful output. It is currently only scripted to run on a single spectrum.}
#' @references {
#' More information and examples of this function can be found in the Soil Spectral Inference training materials and references therein from the Sydney Soil Schools.
#' }
#' @examples 
#' data(soilSpecDat)
#' 
#' #Curation
#' nc <- ncol(soilSpecDat)
#' spectra <- soilSpecDat[, 5:nc] #just get the spectra values
#' wavelength <- seq(350, 2500, by = 1) #wavelength sequence
#' colnames(spectra) <- wavelength #append column names
#' #Trim spectra table to a specified region of interest
#set limits
#' lower<- 2079
#' upper<- 2267
#' specTrim<- trimSpec(spectra, wavlimits=range(lower:upper)) 
#' plot(seq(from=lower, to=upper, by=1),specTrim[1,],type="l") #FIrst Spectrum
#' 
#' #Run the continuum fitting function
#' spec.CR<- chBLC_ext(specTrim[1,], lower, upper) #First spectrum
#' # Draw the continuum
#' lines(seq(from=lower, to=upper, by=1),spec.CR$continuum,col="red") 
#' 
#' #Plot continuum removed spectrum
#' plot(seq(from=lower, to=upper, by=1),spec.CR$c.hull,type="l", xlab="wavelength", ylab="reflectance")
#' polygon(spec.CR$polygon,col='red')
#' @export chBLC_ext
#' @author Mario Fajardo and Brendan Malone

chBLC_ext <- function (spectra, lower = 350, upper = 2500, AUC=TRUE)
{
  interval <- as.numeric(colnames(spectra))
  hull_spectra <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
  tempSpect <-  round(as.matrix(spectra),digits = 4)
  if(any(tempSpect<0)){
    correction<-tempSpect[tempSpect<0][which.min(tempSpect[tempSpect<0]-0)]
    under_0 <- interval[tempSpect<0]
    tempSpect <- tempSpect+abs(correction)
    warning(paste0('Please check your spectra, Wavelengths: ',paste0(under_0,collapse = ' ')," had negative values"))
  } 
  data1 <- sortedXyData(interval, tempSpect)
  data1$y <- round(data1$y,digits = 4)
  # plot(data1$x,data1$y,type='b')
  c_hull <- chull(data1)
  # points(data1$x[c_hull],data1$y[c_hull],col='red',cex=2)
  first_last <- data.frame(x=data1[c(1,nrow(data1)),1],y=data1[c(1,nrow(data1)),2])  
  linear_model <- lm(formula = y~x,data = first_last)
  
  
  if(AUC){ #if AUC select those points higher than the linear approx between extremes
    is_up <- as.logical((data1[c_hull,2]>=round(predict(linear_model,data1[c_hull,]),digits = 4)))  
    c_hull <- c_hull[is_up]
  }
  
  if(!AUC){
    is_down <- as.logical((data1[c_hull,]<=round(predict(linear_model,data1[c_hull,]),digits = 4))[,2])  
    c_hull <- c_hull[is_down]  
  }
  # points(data1[c_hull,],col='blue',pch=16)
  # lines(data1[c_hull,][order(data1[c_hull,1]),],col='blue',pch=16)
  
  linear_approx <- approx(data1[c_hull, ], xout = interval, 
                          method = "linear", ties = "mean")
  linear_approx$y <- round(linear_approx$y,digits = 4)
  if (AUC) {
    hull_spectra[1, ] <- 1-((linear_approx[[2]] - tempSpect)/linear_approx[[2]])
    hull_spectra[, c(1, ncol(hull_spectra))] <- 1.0001 #in order to make the extreme points unique
    hull_data <- data.frame(wave=as.numeric(colnames(tempSpect)),value=t(hull_spectra))
    closing_line <- data.frame(wave=hull_data$wave, value=rep(1.0001, each = length(hull_data$wave)))
  }
  
  if (!AUC) {
    hull_spectra[1, ] <- (tempSpect-linear_approx[[2]])/linear_approx[[2]]
    hull_spectra[, c(1, ncol(hull_spectra))] <- -0.0001 #in order to make the extreme points unique
    hull_data <- data.frame(wave=as.numeric(colnames(tempSpect)),value=t(hull_spectra))
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
  if(AUC){
    
    min_raw <- tempSpect[which(hull_spectra==min(hull_spectra))[1]]  
    min_continuum <- hull_data$value[which(hull_spectra==min(hull_spectra))[1]]
    if(min_raw==0&min_raw==0){
      depth_feature <- 1
    }else{
      depth_feature <- 1-(min_raw/min_continuum)
    }
  }
  
  if(!AUC){
    min_raw <- tempSpect[which(hull_spectra==min(hull_spectra))[1]]
    min_continuum <- hull_data$value[which(hull_spectra==min(hull_spectra))[1]]
    depth_feature <- 1-(min_raw/min_continuum)
    depth_feature <- tempSpect[which(hull_spectra==min(hull_spectra))[1]]/hull_data$value[which(hull_spectra==min(hull_spectra))[1]]
  }
  
  retval <- list(wave = hull_data$wave, 
                 c.hull = hull_spectra, 
                 raw.spec = data1$y, 
                 continuum = linear_approx[[2]], 
                 polygon = hull_polygon,
                 area=area_feature,
                 depth=depth_feature,
                 slope=slope_feature)
  
  return(retval)
}