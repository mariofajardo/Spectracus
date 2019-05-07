#' Train Bootstrapped models
#' 
#'  
#' 
#' @author Mario Fajardo
#' 
#' @param eSpectra A \code{\link{eSpectra}} object.
#' @param bootSample a bootsample obtained with \code{\link{bag_eSpectra}}
#' @param target character vector of variable to be predicted (must be in eSpectra@Properties)
#' 
#' @return A list containing:
#' \item{predicted}{predictions}
#' \item{observed}{observations}
#' \item{acc}{accuracy indices}
#' 
#' 
#' @importFrom stats predict
#' @importFrom pbapply pblapply
#' @importFrom Cubist cubist
#' @importFrom clhs clhs
#' 
#' @export boot_models

boot_models <- function(eSpectra=eSpectra,bootSample=NULL,target=NULL,...){
  
  if(length(bootSample$SamplesIn)==length(eSpectra@ID)) {
  
  #LOOCV Work in progress  
  MODELS <- pblapply(1:ncol(bootSample$BootSamples),function(samples){
    
    validation <- eSpectra@ID[!eSpectra@ID%in%eSpectra[bootSample$BootSamples[,samples]]@ID] 
    model <- cubist(as.data.frame(eSpectra[bootSample$BootSamples[,samples]]@Spectra),eSpectra[bootSample$BootSamples[,samples]]@Properties[,target])
    preds <- predict(model,as.data.frame(eSpectra[validation]@Spectra))
    result <- list(model,preds)
  })
  
  predictions <- rbind(sapply(MODELS,function(x)x[[2]]))
  preds <- as.numeric(predictions)
  obs <- eSpectra[apply(bootSample$BootSamples,2,function(x) bootSample$SamplesIn[!bootSample$SamplesIn%in%x])]@Properties[,target]
    
    
  }else{
    
    
    validation <- eSpectra@ID[!eSpectra@ID%in%bootSample$SamplesIn]
    MODELS <- pbapply(bootSample$BootSamples,2,function(samples){
      model <- cubist(as.data.frame(eSpectra[samples]@Spectra),eSpectra[samples]@Properties[,target])
      preds <- predict(model,as.data.frame(eSpectra[validation]@Spectra))
      result <- list(model,preds)
    })
    
    predictions <- rbind(sapply(MODELS,function(x)x[[2]]))
    preds <- rowMeans(predictions)
    obs <- eSpectra[validation]@Properties[,target]
    
  }
  
  return(list(predicted=predictions,
              observed=obs,
              acc=goof(obs,preds)))
}
