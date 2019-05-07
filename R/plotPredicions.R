#' Plot bootstrapped model predictions
#'  
#' 
#' @author Mario Fajardo
#' 
#' @param models list with three elements (predicted, observed and accuracy data) obtained after using \code{\link{boot_models}} function. 
#' @param formula Should a formula with a linear regression be presented (Logic)
#' @param formulaCoords x and y coordinates of Formula (if formula=TRUE)
#' 
#' @return A \code{\link{ggplot}} object
#' 
#' 
#' @import ggplot2 
#' 
#' @export plotPredicions

plotPredicions <- function(models,formula=FALSE,formulaCoords=NULL){
  
  lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
  }
  
  df=  data.frame(y=rowMeans(models$predicted),x=models$observed)
  
  predictionsData <- data.frame(observations = models$observed,predictions=rowMeans(models$predicted))
  p <- ggplot(predictionsData,aes(observations,predictions))+
    geom_point(size=2,alpha=.8)+
    geom_abline(slope=1)+ 
    geom_smooth(method='lm',formula=y~x,se=FALSE,linetype='dashed')
  p
  
  if (formula) p + geom_text(x = formulaCoords[1], y = formulaCoords[2], label = lm_eqn(df), parse = TRUE,size=5)
  
}
