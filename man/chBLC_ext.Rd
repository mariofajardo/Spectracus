\name{Spectra normalisation via continuum estimation}
\docType{methods}
\alias{chBLC_ext}
\title{Spectra normalisation via continuum estimation}
\description{This function approximates the continuum of a provided spectrum or spectrum region of interest, then derives the normalised spectrum from this continuum i.e. the deviation from estimated continuum.}
\usage{chBLC_ext(spectra, lower=350, upper=2500, specType = 0)}

\arguments{
\item{spectra}{A single spectrum of class \code{data.frame}. Column values represent spectral values per wavelength unit.}
\item{lower}{numeric; The lower wavelength bound of the spectral feature or region of interest.}
\item{upper}{numeric; The upper wavelength bound of the spectral feature or region of interest.}
\item{specType}{A logical that is expressed as either 0 or 1 to indicate whether continuum needs to be fitted to raw data or inverted data.}
 }

\value{Returns a list with 5 elements:
\describe{
\item{wave}{The sequence of spectral wavelengths of the region of interest.}
\item{c.hull}{The normalised or continuum removed spectrum}
\item{raw.spec}{The raw spectrum}
\item{continuum}{The fitted continuum}
\item{polygon}{vertices for constructing the polygon of the normalised spectral feature}
}}


\note{This function works similarly to the other functions in this package where a continuum is fitted to data from which raw spectra values are removed from it. However, this function provides extra useful output. It is currently only scripted to run on a single spectrum.}

\author{Brendan Malone}

\references{
\itemize{
\item More information and examples of this function can be found in the Soil Spectral Inference training materials and references therein from the Sydney Soil Schools.
}
}

\examples{
data(soilSpecDat)

#Curation
nc <- ncol(soilSpecDat)
spectra <- soilSpecDat[, 5:nc] #just get the spectra values
wavelength <- seq(350, 2500, by = 1) #wavelength sequence
colnames(spectra) <- wavelength #append column names

#Trim spectra table to a specified region of interest
#set limits
lower<- 2079
upper<- 2267
specTrim<- trimSpec(spectra, wavlimits=range(lower:upper)) 
plot(seq(from=lower, to=upper, by=1),specTrim[1,],type="l") #FIrst Spectrum

#Run the continuum fitting function
spec.CR<- chBLC_ext(specTrim[1,], lower, upper,specType = 1) #First spectrum
# Draw the continuum
lines(seq(from=lower, to=upper, by=1),spec.CR$continuum,col="red") 

#Plot continuum removed spectrum
plot(seq(from=lower, to=upper, by=1),spec.CR$c.hull,type="l", xlab="wavelength", ylab="reflectance")

}
\keyword{methods}
