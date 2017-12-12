#' Filtering methods for SoilSpectra objects
#' 
#' This method calls different filters on matrices in which each row is a separate spectrum
#' depending on the value of 'type':
#' @rdname  filter_SoilSpectra
#' @param SoilSpectra Object of class \code{\link{SoilSpectra}}.
#' @param type  
#' \itemize{
#' \item 'S-Golay' for Sativsky Golay filter using \code{\link{sgolayfilt}},
#' \item 'Wavelet' for Wavelet decomposition using \code{\link{accessC.wd}},
#' \item 'MSC' for Multiplicative scatter correction, 
#' \item 'SNV' for Standard Normal Variate transformation.
#' \item 'C-hull' for Continuum removal by removing the convex hull
#' }
#' @param If type = 'S-Golay' then n,p,m parameters control \code{\link{sgolayfilt}} window size, polynomial order and derivative order respectively.
#' @param If type = 'Wavelet' res level to be extracted from wavelet decomposition model, see \code{\link{accessC.wd}}
#' @param If type = 'C-hull' specType also needs to be specified, with 0 if the data is in absorbance units or 1 if the data is in reflectance units.
#' @param If type is a character vector with more than one filter type the method will apply all the filters in the order they are specified.
#' @importFrom signal sgolayfilt
#' @importFrom signal sgolayfilt
#' @importFrom wavethresh wd accessC.wd
#' @references Savitzky, Abraham and Golay, M. J. E. (1964) Smoothing and Differentiation of Data by Simplified Least Squares Procedures, Analytical Chemistry, 36, pp. 1627-1629.
#' @author Brendan Malone, Budiman Minasny, Michael Nelson, Sebastian Campbell and Mario Fajardo.
#' @examples 
#' \dontrun{
#' data("SoilSpectraExample")
#' 
#' 
#' plot(SoilSpectraExample@Wavelength,
#'      SoilSpectraExample@Spectra[1,],
#'      main=SoilSpectraExample@Treatments,
#'      ylab=SoilSpectraExample@Type,
#'      xlab=slotNames(SoilSpectraExample)[3])
#' 
#' 
#' #Savitsky-Golay filter
#' SGolay_Spectra <- filter_SoilSpectra(SoilSpectraExample,type = 'S-Golay')
#' 
#' plot(SGolay_Spectra@Wavelength,
#'      SGolay_Spectra@Spectra[1,],
#'      main=SGolay_Spectra@Treatments,
#'      ylab=SGolay_Spectra@Type,
#'      xlab=slotNames(SGolay_Spectra)[3])
#' 
#' #Multiplicative scatter correction
#' 
#' MSC_Spectra <- filter_SoilSpectra(SoilSpectraExample,type = 'MSC')
#' 
#' plot(MSC_Spectra@Wavelength,
#'      MSC_Spectra@Spectra[1,],
#'      main=MSC_Spectra@Treatments,
#'      ylab=MSC_Spectra@Type,
#'      xlab=slotNames(MSC_Spectra)[3])
#' 
#' 
#' 
#' #Standard Normal Variate transform
#' 
#' SNV_Spectra <- filter_SoilSpectra(SoilSpectraExample,type = 'SNV')
#' 
#' plot(SNV_Spectra@Wavelength,
#'      SNV_Spectra@Spectra[1,],
#'      main=SNV_Spectra@Treatments,
#'      ylab=SNV_Spectra@Type,
#'      xlab=slotNames(SNV_Spectra)[3])
#' 
#' 
#' #Convex-hull transform
#' 
#' CHull_Spectra <- filter_SoilSpectra(SoilSpectraExample,type = 'C-hull')
#' 
#' plot(CHull_Spectra@Wavelength,
#'      CHull_Spectra@Spectra[1,],
#'      main=CHull_Spectra@Treatments,
#'      ylab=CHull_Spectra@Type,
#'      xlab=slotNames(CHull_Spectra)[3])
#' }
#' @exportMethod  filter_SoilSpectra


setGeneric("filter_SoilSpectra",
           function(SoilSpectra,type,n=11, p=2, m=0,res=NULL,...)
             {
             standardGeneric('filter_SoilSpectra')
             }
           )

setMethod(f = 'filter_SoilSpectra',
          signature = 'SoilSpectra',
          definition= function(SoilSpectra,type=NULL,n=11, p=2, m=0,res=NULL,...)
            {
    spectra <- SoilSpectra@Spectra
    specType <- SoilSpectra@Type
    
    if(length(type)>1){
      for(treatment in type){
        SoilSpectra <- filter_SoilSpectra(SoilSpectra,treatment)
      }
      return(SoilSpectra)
      
    }else{
      
      if (type =='S-Golay'){
        ## run filter
        if(nrow(spectra)>1){
          sg <- t(apply(spectra, 1, sgolayfilt, n = n, p = p, m = m))
        }else{
          sg <- matrix(sgolayfilt(spectra, n = n, p = p, m = m),ncol = ncol(spectra))
        }  
        
        SoilSpectra@Spectra <- sg
        treatmentDetails <- paste0('S-Golay','n=',n,'p=',p,'m=',m)
        SoilSpectra@Treatments <- c(SoilSpectra@Treatments,treatmentDetails)
        return(SoilSpectra)
      }
      
      if(type=='MSC') {
        #first calculate a mean spectrum.
        meanSpec <- colMeans(spectra)
        mscMat <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
        for (i in 1:nrow(spectra)) {
          # determine the slope and intercept co-efficents
          specLM <- lm(spectra[i, ] ~ meanSpec)
          specCE <- t(as.matrix(specLM$coefficients))
          # Adjust the spectra
          mscMat[i, ] <- t(as.matrix((spectra[i, ] - specCE[1, 1])/specCE[1,2]))
        }
        
        SoilSpectra@Spectra <- mscMat
        treatmentDetails <- 'MSC'
        SoilSpectra@Treatments <- c(SoilSpectra@Treatments,treatmentDetails)
        return(SoilSpectra)
      }
      
      if(type=='SNV'){
        snvMat<-(spectra - rowMeans(spectra))/apply(spectra,1,sd)
        SoilSpectra@Spectra <- snvMat
        treatmentDetails <- 'SNV'
        SoilSpectra@Treatments <- c(SoilSpectra@Treatments,treatmentDetails)
        return(SoilSpectra)
      }
      
      if(type=='Wavelet'){
        nm2<- 2^c(1:100)
        vs<- ncol(spectra)
        if (sum(nm2 == vs) != 1) {
          stop("Error: Number of columns in spectra table needs to equal 2^x")
        }else{
          wave_spectra <- matrix(NA, ncol = 2^res, nrow = nrow(spectra))
          
          for (i in 1:nrow(spectra)){
            wds <- wd(as.matrix(spectra[i, ]), 
                      bc = "symmetric", 
                      filter.number = 10, 
                      family = "DaubExPhase", 
                      min.scale = 2)
            
            wave_spectra[i, ] <- accessC.wd(wds, level = res)
          }
          
          wave_spectra<- as.data.frame(wave_spectra)
          colnames(wave_spectra) <- seq((as.numeric(names(spectra)[1]) + 0.5 * (ncol(spectra)/(2^res))), as.numeric(names(spectra)[length(spectra)]), by = ncol(spectra)/(2^res)) 
        }
        return(wave_spectra)
      }
      
      if(type=='C-hull'){
        c_hull.fix<- function(c_hull){
          cc<-c_hull
          xs<-which(cc == 1)
          ccd<- sort(cc[which(cc == 1):length(cc)])
          if (cc[1]< cc[length(cc)]){
            ccd<-ccd
          }else{
            
            for(f in 1:(xs-1)){
              if(cc[f]< cc[f+1]){
                ccd<- c(ccd,cc[f])
              }else{
                ccd<- c(ccd,cc[f])
                break}
            }
          }
          return(ccd)
        }
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
          if (cc[1]< cc[length(cc)]){
            ccd<-ccd
          }else{
            
            for(f in 1:(xs-1)){
              if(cc[f]< cc[f+1]){
                ccd<- c(ccd,cc[f])
              }else{
                ccd<- c(ccd,cc[f])
                break}
            }
          }
          c_hull <- ccd
          
          ## calculate linear approximation between hull points
          linear_approx <- approx(data1[c_hull,], xout = interval, method = 'linear', ties = 'mean')
          ## calculate the deviation from the convex hull
          if (SoilSpectra@Type=='Absorbance')
          {hull_spectra[i,] <- 1- (( linear_approx[[2]] - tempSpect )/linear_approx[[2]])}
          if (SoilSpectra@Type=='Reflectance')
          {hull_spectra[i,] <- (( linear_approx[[2]] - tempSpect )/linear_approx[[2]])}
        }
        
        SoilSpectra@Spectra <- hull_spectra
        treatmentDetails <- 'C-hull'
        SoilSpectra@Treatments <- c(SoilSpectra@Treatments,treatmentDetails)
        return(SoilSpectra)
      }
        
      }

          }
) 
    


