#' Predict using a fitted bootstrap aggregated PLSR model
#' 
#' #' (NOT YET IMPLEMENTED) 
#'  
#' 
#' @author Budiman Minasny
#' 
#' @param model.bpls bagged PLSR model
#' @param newspec matrix of spectra to be predicted where each column is a wavelength and each row a spectrum
#' @param nbag size of bootstrap sample
#' @param nc number of components
#' 
#' @return A list containing:
#' \item{bag.pred}{predicted values}
#' \item{pred.ave}{mean of predictions}
#' \item{pred.std}{standard deviation of predictions}
#' 
#' @seealso \code{\link{fit_bag_plsr}}
#' 
#' @importFrom stats predict
#' 
#' @export

predict_bag_plsr<-function(model.bpls,newspec,nbag,nc){
  # model.bpls = bagged PLSR model, newspec= spectra to be predicted, nbag, no. of bootstrap, nc= no. of components used in PLSR
  # return variables: bag.pred=bagged predicted values, pred.ave=mean of prediction, pred.std= std. dev of prediction
  # use pls package  
  n<-nrow(newspec)
  pred.v <- matrix(0,nrow=n,ncol=nbag)
  
  for (ibag in 1:nbag) {
    pred.v[,ibag] <- predict(model.bpls[[ibag]], ncomp = nc, newdata = newspec)
  }
  pred.ave<-apply(pred.v, 1, mean)
  pred.std<-apply(pred.v, 1, sd)
  
  return(list(bag.pred=pred.v, pred.ave=pred.ave, pred.std=pred.std)) }