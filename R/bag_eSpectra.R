#' Split eSpectra object into a bootstrapped sample
#' 
#'  
#' 
#' @author Mario Fajardo
#' 
#' @param eSpectra A \code{\link{eSpectra}} object.
#' @param val_split percentage of the dataset that will be used for validation
#' @param mode currently two methods implemented, 'lhc' for latin hypercube sampling and 'random' for random sampling (default). The latin hypercube sample design uses \code{\link{clhs}} function and is applied to the first 5 principal component scores of the eSpectra Spectral information
#' @param lhcIter Latin hypercube iterations (if mode = 'lhc').For testing purposes is recommended to lower this value.  
#' 
#' @return A list containing:
#' \item{BootSamples}{Bootstrapped samples IDs}
#' \item{SamplesIn}{Sample IDs used for training}
#' 
#' 
#' @importFrom stats predict
#' @importFrom clhs clhs
#' 
#' 
#' @examples 
#' \dontrun{
#' data("eSpectraExample")
#'
#' Bag_Sample <- bag_eSpectra(eSpectra=eSpectraExample,val_split = 0.25,mode = 'random',iters = 10)
#' 
#' Bag_Sample
#' 
#' }
#' 
#' @export bag_eSpectra

bag_eSpectra <- function(eSpectra= NULL,val_split=NULL,mode='random',iters=50,lhcIter=15000,...){
  
  # if(validation=='LOOCV'){
  # 
  #   if (mode=='random'){
  #     bootSamples <- sapply(1:iters,function(x){
  #       set.seed(x)
  #       train <- sample(eSpectra@ID,size=length(eSpectra@ID)-1)
  #     })
  #     return(list(BootSamples=bootSamples,SamplesIn=eSpectra@ID))
  #   }
  # 
  #   if(mode=='lhc'){
  #     comps <- as.data.frame(prcomp(eSpectra@Spectra,center = T,scale. = T)$x[,1:5])
  #     bootSamples <- sapply(1:iters,function(x){
  #       set.seed(x)
  #       train <- clhs(comps,size=length(eSpectra@ID)-1,iter=lhcIter)
  #       eSpectra@ID[train]
  #     })
  #     return(list(BootSamples=bootSamples,SamplesIn=eSpectra@ID))
  #   }
  # }
  
    if (mode=='random'){
      
      SamplesIn <- sample(eSpectra@ID,size=length(eSpectra@ID)*(1-val_split))
      bootSamples <- sapply(1:iters,function(x){
        set.seed(x)
        train <- sample(SamplesIn,size=length(SamplesIn)-(length(SamplesIn)*.01))
      })
      return(list(BootSamples=bootSamples,SamplesIn=SamplesIn))
    }
    
    if(mode=='lhc'){
      
      
      SamplesIn <- sample(eSpectra@ID,size=length(eSpectra@ID)*(1-val_split))
      val <- eSpectra@ID[!eSpectra@ID%in%SamplesIn]
      comps <- as.data.frame(prcomp(eSpectra@Spectra[eSpectra@ID%in%SamplesIn,],center = T,scale. = T)$x[,1:5])
      bootSamples <- sapply(1:iters,function(x){
        set.seed(x)
        samplesInforTrain <- sample(SamplesIn,size=length(SamplesIn)-(length(SamplesIn)*.1))
        train <- clhs(comps,size=length(samplesInforTrain),iter=lhcIter)
        eSpectra@ID[train]
      })
      
      return(list(BootSamples=bootSamples,SamplesIn=SamplesIn))
      
    }    
    
    

}