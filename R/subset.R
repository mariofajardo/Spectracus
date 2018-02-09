#' SoilSpectra subset method


setMethod(f = "[", 
          signature=c("SoilSpectra"),
          definition=function(x,i,j,...,drop=T)
          {
            if(missing(i)) i <- as.character(x@ID)
            if(missing(j)) j <- as.character(x@Bands)
            
            if(is.numeric(i)) logic_i <- 1:length(x@ID)%in%i
            if(is.character(i)) logic_i <- x@ID%in%i
            if(is.logical(i))   logic_i <- i
            
            if(is.numeric(j)) logic_j <- x@Bands%in%x@Bands[j]
            if(is.character(j)) logic_j <- x@Bands%in%j   
            if(is.logical(j))   logic_j <- j
            
            initialize(x,
                      Instrument=x@Instrument,
                      Spectra=if(length(i)>1){
                        
                        if(length(j)==1) {
                          t(matrix(x@Spectra[logic_i,logic_j]))
                        }else{
                          (x@Spectra[logic_i,logic_j])
                        }
                        
                        }else{
                          
                          t(matrix(x@Spectra[logic_i,logic_j]))
                          
                          },
                      
                      Bands=x@Bands[logic_j],
                      RowsAreSpectra=x@RowsAreSpectra,
                      Type=x@Type,
                      ID=x@ID[logic_i],
                      Properties=x@Properties[logic_i,],
                      Treatments=x@Treatments)
          })
