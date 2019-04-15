#' Create Munsell treemap (NOT YET IMPLEMENTED)
#' 
#' Creates a treemap with frequency of munsell colours as the size of the tree rectangles using \code{treemap} and actual Munsell colours as the colours of rectangles using \code{ggplot2}.
#' 
#' @author Michael Nelson, Sebastian Campbell, Mario Fajardo
#' 
#' @importFrom munsell rgb2mnsl
#' @importFrom munsell mnsl
#' @importFrom plyr splat
#' @importFrom plyr adply
#' @importFrom plyr llply
#' @importFrom treemap treemap 
#' @import ggplot2 

#' 
#' @param spectra dataframe or matrix where each row is a spectrum and each column is a wavelength
#' @param wavelengths integer wavelengths corresponding to the columns of \code{spectra}
#' @param coltext logical, whether or not to plot Munsell colours as text
#' @param numtext logical, whether or not to plot colour frequencies
#' @param otherArgs list of additional elements to be added to ggplot2 call
#' @param textrange vector of length two indicating minimum and maximum text size respectively
#' @export munsell_treemap

munsell_treemap <- function(spectra, wavelengths, coltext=TRUE, numtext=FALSE, otherArgs=NULL, textrange=c(3, 10)){
  .e<-environment()
  treemap_data <- munsell_tm(spectra, wavelengths)
  treemap_plot <- ggplot(treemap_data, aes(xmin=x0, xmax=x0+w, ymin=y0, ymax=y0+h), environment=.e)+
    geom_rect(aes(fill=hex))+
    scale_fill_identity()+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())
    
  if(coltext){
    upline <- if(numtext) 0.66 else 0.5    
    treemap_plot <- treemap_plot + 
    geom_text(aes(label=munsell, x=x0+0.5*w, y=y0+upline*h, size=log10(Freq)), show_guide=FALSE)+
    scale_size_continuous(range=textrange)
  }
  if(numtext){
    lowline <- if(coltext) 0.33 else 0.5
    treemap_plot <- treemap_plot +
    geom_text(aes(label=Freq, x=x0+0.5*w, y=y0+lowline*h, size=log10(Freq)), 
                show_guide=FALSE)
  }
    
  print(treemap_plot+otherArgs)
}


spectra_to_RGB <- function(.spectra, all_wavelengths,rgb_list = list(blue = list(.interval=450:520), red = list(.interval=600:690), green=list(.interval=520:600)), ...) {
  in_interval <- function(.all, .interval,...){
    ## index for an interval
    ## a function to subset a particular waveband interval
    .in_interval = .all %in% .interval
  }
  mean_interval <- function(.data, .index){
    ## returns the mean for given indices
    mean(.data[.index])
  }
  
  interval_list <- llply(rgb_list, splat(in_interval), .all=all_wavelengths)
  ##
  ## get the average in these subsets
  rgb_values <- lapply(interval_list, mean_interval, .data=.spectra)
  ##
  ## convert to colour
  colour <- with(rgb_values, rgb(red, green, blue))
  ##
  ## return data frame
  with(rgb_values, data.frame(red, green, blue, colour))
}
munsell_tm <- function(spectra, wavelengths){
  
  spectra <- as.matrix(spectra)
  
  rgb_colours <- adply(spectra,1, spectra_to_RGB, all_wavelengths = wavelengths)
  munsell <- splat(function(red,green,blue, ...){rgb2mnsl(R=red, G=green, B=blue)})(rgb_colours)
  
  munsell_table <- as.data.frame(table(munsell))
  
  pdf(file=NULL)
  raw_tmdata <- treemap::treemap(munsell_table, index="munsell", vSize="Freq")[[1]]
  dev.off()
  
  ordered_munsell <- sapply(raw_tmdata[,-c(ncol(raw_tmdata)-0:4)], levels)
  names(ordered_munsell) <- NULL
  
  rect_coords <- data.frame(munsell=ordered_munsell[[1]], 
                            raw_tmdata[,c(ncol(raw_tmdata)-4:0)])
  rect_coords$hex <- munsell::mnsl(rect_coords$munsell)
  rect_coords <- merge(rect_coords, munsell_table)
  
  rect_coords
}
