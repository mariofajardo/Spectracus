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
filter_spectra <-function(spectra,type,n=11, p=2, m=0,res=NULL){
  
  if (type =='S-Golay'){
      spectra <- as.matrix(spectra)
    ## run filter
    sg <- aaply(spectra, 1, sgolayfilt, n = n, p = p, m = m)
    ##
    ## arrange appropriately if a single sample
    if(nrow(spectra) == 1){sg <- matrix(sg, dim(spectra))}
    ## return data frame
    colnames(sg) <- colnames(spectra)
    return(sg)
  } 

  if(type=='MSC') {
      spectra<-as.matrix(spectra)
      #calculate the mean spectrum
      meanSpec<-as.matrix(colMeans(spectra))
      mscMat<- matrix(NA,ncol=ncol(spectra),nrow=nrow(spectra))
      for (i in 1:nrow(spectra)){
        #determine the slope and intercept co-efficents 
        specLM<- lm(spectra[i,]~ meanSpec)
        specCE<-t(as.matrix(specLM$coefficients))
        #Adjust the spectra
        mscMat[i,]<- t(as.matrix((spectra[i,]-specCE[1,1])/specCE[1,2]))}
      colnames(mscMat) <- colnames(spectra)
      return(mscMat)
    
  }
  
  if(type=='SNV'){
    
      spectra<-as.matrix(spectra)
      snvMat<-(spectra - rowMeans(spectra))/apply(spectra,1,sd)
      attributes(snvMat) <- attributes(spectra)
      return(snvMat)
  }

  if(type=='Wavelet'){
      spectra<-as.matrix(spectra)
      wave_spectra<- matrix(NA,ncol=2^res,nrow=nrow(spectra))
      for (i in 1:nrow(spectra)){
        wds<-wd(as.matrix(spectra[i,]),bc="symmetric",filter.number = 10, family = 'DaubExPhase', min.scale = 2)
        wave_spectra[i,]<- accessC.wd(wds, level=res)}
      colnames(wave_spectra) <- seq((404 + 0.5*(2048/(2^res))),2451, by=2048/(2^res))
      return(wave_spectra)
  }
  
  if(type=='C-hull'){
      spectra<-as.matrix(spectra)
      interval <- seq_len(ncol(spectra))
      hull_spectra <- matrix(NA,ncol=ncol(spectra),nrow=nrow(spectra))
      for (i in seq_len(nrow(spectra))){
        tempSpect <- as.matrix(spectra[i,])
        data1 <- sortedXyData(interval, tempSpect)
        ## calculate convex hull
        c_hull <- chull(data1)
        ## get the appropriate region: the points of the polygon over the spectra
        
        # Create vector which wraps around
        c_hull <- c(c_hull, c_hull)
        # remove all points before the first one.
        c_hull <- c_hull[which.min(c_hull):length(c_hull)]
        # Go until the first end
        c_hull <- c_hull[1:which.max(c_hull)]
        
        ## calculate linear approximation between hull points
        linear_approx <- approx(data1[c_hull,], xout = interval, method = 'linear', ties = 'mean')
        ## calculate the deviation from the convex hull
        hull_spectra[i,] <- ( linear_approx[[2]] - tempSpect )/linear_approx[[2]]}
      colnames(hull_spectra) <- colnames(spectra)
      return(hull_spectra)
      }

}

