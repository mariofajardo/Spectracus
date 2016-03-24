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
#' @param specType A logical that is expressed as either 0 or 1 to indicate whether continuum needs to be fitted to raw data or inverted data.}
#' 
#' @return Returns a list with 5 elements:
#' \itemize{
#' \item{wave}{ The sequence of spectral wavelengths of the region of interest.}
#' \item{c.hull}{ The normalised or continuum removed spectrum}
#' \item{raw.spec}{ The raw spectrum}
#' \item{continuum}{ The fitted continuum}
#' \item{polygon}{ Vertices for constructing the polygon of the normalised spectral feature}
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
#' spec.CR<- chBLC_ext(specTrim[1,], lower, upper,specType = 1) #First spectrum
#' # Draw the continuum
#' lines(seq(from=lower, to=upper, by=1),spec.CR$continuum,col="red") 
#' 
#' #Plot continuum removed spectrum
#' plot(seq(from=lower, to=upper, by=1),spec.CR$c.hull,type="l", xlab="wavelength", ylab="reflectance")
#' @export chBLC_ext
#' @author Brendan Malone

chBLC_ext<- function(spectra, lower=350, upper=2500, specType = 0){
  interval<- c(1:ncol(spectra))
  hull_spectra<- matrix(NA,ncol=ncol(spectra),nrow=nrow(spectra))
  tempSpect= as.matrix(spectra)
  data1 <- sortedXyData(interval, tempSpect)
  
  ## calculate convex hull
  c_hull <- chull(data1)
  ## get the appropriate region
  cc<-c_hull
  xs<-which(cc == 1)
  ccd<- cc[which(cc == 1):length(cc)]
  ccd.m<- matrix(NA, nrow=1, ncol=length(ccd))
  for(f in 1:(length(ccd)-1)){
    if(ccd[f] < ccd[f+1]){ccd.m[1,f]<- ccd[f]}
    else {ccd.m[1,f]<- ccd[f]
    break}}
  if(f == (length(ccd)-1)){ccd.m[1,ncol(ccd.m)]<- ccd[length(ccd)]}
  tempX<- c(ccd.m)[!is.na(c(ccd.m))]
  #Left side
  if(c_hull[1] > c_hull[length(c_hull)]){lss<- cc[1:(xs-1)]
  slss<- lss[which(lss >= lss[1])]
  c_hull<- sort(c(tempX,slss))} else {c_hull<-tempX}
  
  ## calculate linear approximation between hull points
  linear_approx <- approx(data1[c_hull,], xout = interval, method = 'linear', ties = 'mean')
  ## calculate the deviation from the convex hull
  if (specType == 1)
  {hull_spectra[1,] <- 1- (( linear_approx[[2]] - tempSpect )/linear_approx[[2]])}
  if (specType == 0)
  {hull_spectra[1,] <- (( linear_approx[[2]] - tempSpect )/linear_approx[[2]])}
  colnames(hull_spectra) <- colnames(spectra)
  # construct polygon
  aaa<- hull_spectra
  aaa[,c(1,ncol(aaa))]<- 1.0001
  bbb<- cbind(seq(lower,upper,by =1),t(aaa))
  ccc<- cbind(seq(lower,upper,by =1),rep(c(1.0001),each=length(seq(lower,upper,by=1))))
  ddd<- as.matrix(rbind(bbb,ccc))
  retval<- list(wave=seq(lower,upper,by=1),c.hull=hull_spectra, raw.spec=data1$y,continuum=linear_approx[[2]] ,polygon=ddd)
  return(retval)}