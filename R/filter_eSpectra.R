#' Filtering methods for eSpectra objects
#' 
#' This method calls different filters on matrices in which each row is a separate spectrum
#' depending on the value of 'type':
#' @rdname  filter_eSpectra
#' @param eSpectra Object of class \code{\link{eSpectra}}.
#' @param type  
#' \itemize{
#' \item 'S-Golay' for Sativsky Golay filter using \code{\link{sgolayfilt}},
#' \item 'Wavelet' for Wavelet decomposition using \code{\link{accessC.wd}},
#' \item 'MSC' for Multiplicative scatter correction, 
#' \item 'SNV' for Standard Normal Variate transformation.
#' \item 'C-hull' for Continuum removal by removing the convex hull
#' \item 'CompSpec' for reducing dimensions of a eSpectra by averaging values inside a window of a specified size
#' }
#' @param If type = 'S-Golay' then n,p and m parameters control \code{\link{sgolayfilt}} window size, polynomial order and derivative order respectively.
#' @param If type = 'Wavelet' res level to be extracted from wavelet decomposition model, see \code{\link{accessC.wd}}
#' @param If type = 'C-hull' specType also needs to be specified, with 0 if the data is in absorbance units or 1 if the data is in reflectance units.
#' @param If type = 'CompSpec' then window parameter will be the size over which the spectra will be averaged
#' @param If type is a character vector with more than one filter type the method will apply all the filters in the order they are specified. Extra parameters will be passed to ...=
#' @importFrom signal sgolayfilt
#' @importFrom signal sgolayfilt
#' @importFrom prospectr spliceCorrection
#' @importFrom wavethresh wd accessC.wd
#' @references Savitzky, Abraham and Golay, M. J. E. (1964) Smoothing and Differentiation of Data by Simplified Least Squares Procedures, Analytical Chemistry, 36, pp. 1627-1629.
#' @author Mario Fajardo, Brendan Malone, Budiman Minasny, Michael Nelson, Sebastian Campbell.
#' @examples 
#' \dontrun{
## Not run: 
#' data("eSpectraExample")
#' 
#' 
#' plot(eSpectraExample[1,as.character(350:2500)]) 

#' #Savitsky-Golay filter
#' SGolay_Spectra <- filter_eSpectra(eSpectraExample,type = 'S-Golay')
#' 

#' plot(SGolay_Spectra[1,as.character(350:2500)])
#' 
#' #Wavelet smoothing
#' 
#' plot(eSpectraExample[1,as.character(350:2397)],ylim=c(0,3))
#' 
#' points(filter_eSpectra(eSpectraExample[1,as.character(350:2397)],type = 'Wavelet',res=9),col='red')
#' points(filter_eSpectra(eSpectraExample[1,as.character(350:2397)],type = 'Wavelet',res=8),col='blue')
#' points(filter_eSpectra(eSpectraExample[1,as.character(350:2397)],type = 'Wavelet',res=7),col='green')
#' 
#' #Multiplicative scatter correction
#' 
#' MSC_Spectra <- filter_eSpectra(eSpectraExample,type = 'MSC')
#' 
#' plot(MSC_Spectra[1,as.character(350:2500)])
#'  
#' #Standard Normal Variate transform
#' 
#' SNV_Spectra <- filter_eSpectra(eSpectraExample,type = 'SNV')
#' 
#' plot(SNV_Spectra[1,as.character(350:2500)])
#' 
#' 
#' #Convex-hull transform
#' 
#' CHull_Spectra <- filter_eSpectra(eSpectraExample,type = 'C-hull')
#' 
#' plot(CHull_Spectra[1,as.character(350:2500)])
#' 
#' 
#' #CompSpec reduction
#' 
#' CompSpec_Spectra <- filter_eSpectra(eSpectraExample,type = 'CompSpec',window=9)
#' 
#' plot(eSpectraExample[1,as.character(350:2500)])
#' points(CompSpec_Spectra[1,as.character(350:2500)],col='red',lwd=2)
#' }
#' 
#' @exportMethod  filter_eSpectra


setGeneric("filter_eSpectra",
           function(eSpectra,type,n=11, p=2, m=0,window=NULL,res=NULL,X = NULL,wav = NULL,splice = NULL,interpol.bands = NULL)
             {
             standardGeneric('filter_eSpectra')
             }
           )

setMethod(f = 'filter_eSpectra',
          signature = 'eSpectra',
          definition= function(eSpectra,type=NULL,n=11, p=2, m=0,window=1,res=NULL,X,wav,splice=c(1000,1830),interpol.bands=10)
            {
            spectra <- eSpectra@Spectra
            specType <- eSpectra@Type
            
            if(length(type)>1){
              for(treatment in type){
                eSpectra <- filter_eSpectra(eSpectra,treatment,n=n, p=p, m=m,window=window,X = X,wav = wav,splice = splice,interpol.bands = interpol.bands)
                }
              return(eSpectra)
              }else{
                
                if (type =='S-Golay')
                  {
        ## run filter
        if(nrow(spectra)>1){
          sg <- t(apply(spectra, 1, sgolayfilt, n = n, p = p, m = m))
        }else{
          sg <- matrix(sgolayfilt(spectra, n = n, p = p, m = m),ncol = ncol(spectra))
        }  
        
        eSpectra@Spectra <- sg
        treatmentDetails <- paste0('S-Golay ','n=',n,'p=',p,'m=',m)
        eSpectra@Treatments <- c(eSpectra@Treatments,treatmentDetails)
        return(eSpectra)
      }
                
                if(type=='MSC') 
                  {
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
        
        eSpectra@Spectra <- mscMat
        treatmentDetails <- 'MSC'
        eSpectra@Treatments <- c(eSpectra@Treatments,treatmentDetails)
        return(eSpectra)
      }
                
                if(type=='SNV')
                  {
        snvMat<-(spectra - rowMeans(spectra))/apply(spectra,1,sd)
        eSpectra@Spectra <- snvMat
        treatmentDetails <- 'SNV'
        eSpectra@Treatments <- c(eSpectra@Treatments,treatmentDetails)
        return(eSpectra)
      }
                
                if(type=='Wavelet')
                  {
        nm2<- 2^c(1:100)
        vs<- ncol(spectra)
        if (sum(nm2 == vs) != 1) {
          stop("Error: Number of columns in spectra table needs to be power of two e.g., 32,64,128,512,1024,2048, etc.")
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
          
        }
        
        eSpectra@Spectra <- wave_spectra
        treatmentDetails <- paste0('Wavelet',' res=', res)
        eSpectra@Treatments <- c(eSpectra@Treatments,treatmentDetails)
        eSpectra@Bands <- as.character(seq((as.numeric(eSpectra@Bands[1]) + 0.5 * (length(eSpectra@Bands)/(2^res))), as.numeric(rev(eSpectra@Bands)[1]), by = length(eSpectra@Bands)/(2^res)))
        
        
        return(eSpectra)
      }
                
                if(type=='C-hull')
                  {
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
          if (eSpectra@Type=='Absorbance')
          {hull_spectra[i,] <- 1- (( linear_approx[[2]] - tempSpect )/linear_approx[[2]])}
          if (eSpectra@Type=='Reflectance')
          {hull_spectra[i,] <- (( linear_approx[[2]] - tempSpect )/linear_approx[[2]])}
        }
        
        eSpectra@Spectra <- hull_spectra
        treatmentDetails <- 'C-hull'
        eSpectra@Treatments <- c(eSpectra@Treatments,treatmentDetails)
        return(eSpectra)
      }
                
                if(type=='CompSpec')
                  {
        
          if(ncol(spectra)%%window != 0){stop("Error: Pick a more compatable window size, the window needs to be a divisor of ",length(eSpectra@Bands), "(total number of bands)")
            
            } else {
              compMat <- matrix(NA, ncol = (ncol(spectra))/window, nrow = nrow(spectra))
              cc <- 1
              for (i in 1:ncol(compMat)) {
                
                if(nrow(compMat)==1) compMat[, i] <- mean(spectra[, cc:(cc + (window - 1))])
                else compMat[, i] <- rowMeans(spectra[, cc:(cc + (window - 1))])
                cc <- cc + window}
                colab = seq(as.numeric(eSpectra@Bands[1]),as.numeric(tail(eSpectra@Bands)[6]), by = window)
                }
          
          eSpectra@Spectra <- compMat
          eSpectra@Bands <- as.character(colab)
          treatmentDetails <- paste0('CompSpec',' window=',window)
          eSpectra@Treatments <- c(eSpectra@Treatments,treatmentDetails)
          return(eSpectra)
                }
                
                if(type=='toAbsorbance')
                {
                  if(eSpectra@Type=='Reflectance'){
                    eSpectra@Spectra <- log(1/eSpectra@Spectra)
                    eSpectra@Type <- 'Absorbance'
                    return(eSpectra)

                  }else{
                    return(eSpectra)  

                    
                }}
                
                if(type=='toReflectance')
                {
                  if(eSpectra@Type=='Absorbance'){
                    eSpectra@Spectra <- log(1/eSpectra@Spectra)
                    eSpectra@Type <- 'Reflectance'
                    return(eSpectra)

                  }else{
                    return(eSpectra)  

                  }
                  
                }
                  
                if(type=='Splice')
                {
                  
                  eSpectra@Spectra <- spliceCorrection(X = eSpectra@Spectra,wav = as.numeric(eSpectra@Bands),splice = splice,interpol.bands = interpol.bands)
                  treatmentDetails <- 'SpliceCorrection'
                  eSpectra@Treatments <- c(eSpectra@Treatments,treatmentDetails)
                  return(eSpectra)
                }
                
              }
            
            }
      ) 
    


