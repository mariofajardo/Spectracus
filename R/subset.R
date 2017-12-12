#' SoilSpectra subset method


setMethod(f = "[", 
          signature=c("SoilSpectra"),
          definition=function(x,i,j,...,drop=TRUE)
          {
           initialize(x,
                      Instrument=x@Instrument,
                      Spectra=if(length(i)>1) {x@Spectra[i,x@Wavelength%in%j]}else{t(matrix(x@Spectra[i,x@Wavelength%in%j]))},
                      Wavelength=x@Wavelength[x@Wavelength%in%j],
                      Range=x@Range,
                      Wavenumber=x@Wavenumber[x@Wavelength%in%j],
                      RowsAreSpectra=x@RowsAreSpectra,
                      Type=x@Type,
                      ID=x@ID[i],
                      Properties=x@Properties[i,],
                      Treatments=x@Treatments)
          })
