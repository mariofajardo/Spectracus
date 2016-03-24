#Generic Function for:
#1. Fitting continuum from spectral hull features
#2. Outputting continuum removed or normalised spectra
#3. Estimating normalised hull feature polygon vertices


# general notes:
# This function is useful in studies for example where spectral features of phenomena need to be investigated
# Selecting the convex hull features in the right order has been an ongoing teething issue.


#Usage:
#Spectra: the spectral data that needs to be processed
#Lower: lower bound of the wavlength range
#Upper: Upper bound of the wavelength range
#SpecType: Rule of thumb is 1 for reflectance data, 0 for absorbance data or otherwise.

#Output
#Returns a list with 5 elements
#wave: The sequence of spectral wavelengths of the region of interest
#c.hull: The normalised or continuum removed spectra
#raw.spec: The raw spectra
#continuum: The fitted continuum
#polygon: vertices for constructing the polygon of the normalised spectral feature

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