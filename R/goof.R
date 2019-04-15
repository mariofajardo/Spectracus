#' Goodness of fit 
#' 
#' #' (NOT YET IMPLEMENTED) 
#'  
#' 
#' @author Brendan Malone
#' 
#' @param observed numeric vector of observed values
#' @param predicted numeric vector of predicted values
#' @param coefficient character vector with performance coefficients as follows:
#' \itemize{
#' \item{'R2'}{ for coefficient of determination}
#' \item{'concordance'}{ for coefficient of concordance as described in Lin, (1989)}
#' \item{'MSE'}{ for mean squared error}
#' \item{'RMSE'}{ for root mean squared error}
#' \item{'bias'}{ for bias coefficient \emph{i.e., bias = mean predicted values - mean observed values}}
#' \item{'MSEc'}{ for bias corrected mean squared error}
#' \item{'RMSEc'}{ for bias corrected root mean squared error}
#' \item{'RPD'}{ for ratio of performance to deviation}
#' \item{'RPIQ'}{ for ratio of performance to interquartile distance}
#' }
#' @param ... arguments passed to \code{\link{plot}}
#' @references {Lin, L., 1989. A Concordance Correlation Coefficient to Evaluate Reproducibility. Biometrics 45(1), 255-268.}
#' @examples
#' observed <- c(1,5,6,8,9,20)
#' predicted <- c(2,5,6,7,9,22)
#' goof(observed,predicted)
#' goof(observed,predicted,coefficient = c('R2','RMSE','RPD'),cex=3,pch=19,col='blue')
#' @export

goof <- function(observed,
                 predicted,
                 coefficient=c('R2','concordance','MSE','RMSE','bias','MSEc','RMSEc','RPD','RPIQ'),
                 plot=TRUE,...){
  
  if(any(!coefficient%in%c('R2','concordance','MSE','RMSE','bias','MSEc','RMSEc','RPD','RPIQ'))) stop('Please choose a valid coefficient')
  
  # Coefficient of determination
  rLM <- lm(predicted ~ observed)
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  
  # Standard error of prediction ^2
  SEP2 <- mean((observed - predicted)^2)
  
  # Standard error of prediction
  SEP <- sqrt(SEP2)
  
  #Bias
  bias <- mean(predicted) - mean(observed)
  
  # residual  variance
  SEP2c <- sum(((predicted - bias - observed)^2) / length(observed))
  SEPc <- sqrt(SEP2c)
  
  # ratio of performance to deviation
  RPD <- sd(observed) / SEP
  
  # Ratio of performance to interquartile distance
  IQ <- c(quantile(observed))[3] - c(quantile(observed))[2]
  RPIQ <- IQ / SEP
  
  # Concordance
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed-mx) * (predicted-my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  
  if(plot){
  plot(observed, 
       predicted,
       ylim=c(min(c(observed,predicted)),max(c(observed,predicted))),
       xlim = c(min(c(observed,predicted)),max(c(observed,predicted))),
       asp=1,
       ...)
    
  abline(a = 0, b = 1, col = "brown4")
  }
  coefs_tmp <- data.frame(R2=R2, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, 
                          MSEc=SEP2c,RMSEc=SEPc, RPD=RPD, RPIQ=RPIQ, row.names=NULL)
  

  gf <- data.frame(coefs_tmp[,coefficient])
  
  gf
}
