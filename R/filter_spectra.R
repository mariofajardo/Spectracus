#' Filtering functions for spectra
#' 
#' This function calls different filters on matrices in which each row is a separate spectrum
#' depending on the value of 'type':
#' \itemize{
#' \item 'S-Golay' for Sativsky Golay filter using \code{\link{sgolayfilt}},
#' \item 'Wavelet' for Wavelet decomposition using \code{\link{accessC.wd}},
#' \item 'MSC' for Multiplicative scatter correction, 
#' \item 'SNV' for Standard Normal Variate transformation.
#' \item 'C-hull' for Continuum removal by removing the convex hull
#' }
#' 
#' @param spectra matrix where each row is a spectrum and each column is a separate wavelength. It is important that these wavelengths be at least very nearly evenly spaced.
#' @param If type = 'S-Golay' then n,p,m parameters control \code{\link{sgolayfilt}} window size, polynomial order and derivative order respectively.
#' @param If type = 'Wavelet' res level to be extracted from wavelet decomposition model, see \code{\link{accessC.wd}}
#' @importFrom plyr aaply
#' @importFrom signal sgolayfilt
#' @importFrom wavethresh wd accessC.wd
#' @export filter_spectra
#' @references Savitzky, Abraham and Golay, M. J. E. (1964) Smoothing and Differentiation of Data by Simplified Least Squares Procedures, Analytical Chemistry, 36, pp. 1627-1629.
#' @author Brendan Malone, Budiman Minasny, Michael Nelson, Sebastian Campbell and Mario Fajardo.
filter_spectra <-function(spectra,type,n=11, p=2, m=0,res=NULL, specType=NULL){
  
  if (type =='S-Golay'){
    spectra <- as.matrix(spectra)
    ## run filter
    sg <- aaply(spectra, 1, sgolayfilt, n = n, p = p, m = m)
    ## arrange appropriately if a single sample
    if (nrow(spectra) == 1) {
      sg <- matrix(sg, dim(spectra))}
    ## return data frame
    sg<- as.data.frame(sg)
    colnames(sg) <- colnames(spectra)
    return(sg)
  } 
  
  if(type=='MSC') {
    #first calculate a mean spectrum.
    meanSpec <- as.matrix(colMeans(spectra))
    mscMat <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
    spectra <- as.matrix(spectra)
    for (i in 1:nrow(spectra)) {
      # determine the slope and intercept co-efficents
      specLM <- lm(spectra[i, ] ~ meanSpec)
      specCE <- t(as.matrix(specLM$coefficients))
      # Adjust the spectra
      mscMat[i, ] <- t(as.matrix((spectra[i, ] - specCE[1, 1])/specCE[1,2]))}
    mscMat<- as.data.frame(mscMat)
    colnames(mscMat) <- colnames(spectra)
    return(mscMat)}
  
  if(type=='SNV'){
    spectra <- as.matrix(spectra)
    snvMat <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
    for (i in 1:nrow(spectra)) {
      snvMat[i, ] <- (spectra[i, ] - mean(spectra[i, ]))/sd(spectra[i,])}
    snvMat<- as.data.frame(snvMat)
    colnames(snvMat) <- colnames(spectra)
    return(snvMat)
  }
  
  if(type=='Wavelet'){
    nm2<- 2^c(1:100)
    vs<- ncol(spectra)
    if (sum(nm2 == vs) != 1) {stop("Error: Number of columns in spectra table needs to equal 2^x")} else {wave_spectra <- matrix(NA, ncol = 2^res, nrow = nrow(spectra)) 
    for (i in 1:nrow(spectra)){ 
      wds <- wd(as.matrix(spectra[i, ]), bc = "symmetric", filter.number = 10, family = "DaubExPhase", min.scale = 2) 
      wave_spectra[i, ] <- accessC.wd(wds, level = res)} 
    wave_spectra<- as.data.frame(wave_spectra)
    colnames(wave_spectra) <- seq((as.numeric(names(spectra)[1]) + 0.5 * (ncol(spectra)/(2^res))), as.numeric(names(spectra)[length(spectra)]), by = ncol(spectra)/(2^res)) }
    return(wave_spectra)
  }
  
  if(type=='C-hull'){
    c_hull.fix<- function(c_hull){
      cc<-c_hull
      xs<-which(cc == 1)
      ccd<- sort(cc[which(cc == 1):length(cc)])
      if (cc[1]< cc[length(cc)])
      {ccd<-ccd} else
      {for(f in 1:(xs-1)){
        if(cc[f]< cc[f+1])
        {ccd<- c(ccd,cc[f])}
        else {ccd<- c(ccd,cc[f])
        break}}}
      return(ccd)}
    
    interval<- c(1:ncol(spectra))
    hull_spectra<- matrix(NA,ncol=ncol(spectra),nrow=nrow(spectra))
    for (i in 1:nrow(spectra)){
      tempSpect= as.matrix(spectra[i,])
      data1 <- sortedXyData(interval, tempSpect)
      
      ## calculate convex hull
      c_hull <- chull(data1)
      
      ## get the appropriate region
      cc<-c_hull
      xs<-which(cc == 1)
      ccd<- sort(cc[which(cc == 1):length(cc)])
      if (cc[1]< cc[length(cc)])
      {ccd<-ccd} else
      {for(f in 1:(xs-1)){
        if(cc[f]< cc[f+1])
        {ccd<- c(ccd,cc[f])}
        else {ccd<- c(ccd,cc[f])
        break}}}
      c_hull <- ccd
      
      ## calculate linear approximation between hull points
      linear_approx <- approx(data1[c_hull,], xout = interval, method = 'linear', ties = 'mean')
      ## calculate the deviation from the convex hull
      if (specType == 1)
      {hull_spectra[i,] <- 1- (( linear_approx[[2]] - tempSpect )/linear_approx[[2]])}
      if (specType == 0)
      {hull_spectra[i,] <- (( linear_approx[[2]] - tempSpect )/linear_approx[[2]])}
    }
    hull_spectra<- as.data.frame(hull_spectra)
    colnames(hull_spectra) <- colnames(spectra)
    return(hull_spectra)
  }
  
}
