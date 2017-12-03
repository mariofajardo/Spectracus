\name{Soil and associated Vis-NIR spectral data}
\docType{data}
\alias{soilSpecDat}
\title{Common data set for soil spectroscopy studies}
\description{A \code{dataframe} of measured soil properties with companion measured and unprocessed vis-nir spectra. The number of soil and spectral observations amounts to 125.}
\usage{data(soilSpecDat)}

\format{
\code{soilSpecDat} is 125 row \code{dataframe} with the first 4 columns indicating observation identifier, measured clay content, total carbon and pH. The remaining columns are visible and near infrared reflectence data. The spectral range is 350-2500nm with 1nm resolution. }

\details{This is a common soil and associated spectral data set that would be used for soil spectroscopy studies}

\note{This data set can be used as a companion to functions contained in the \code{\link{spectroscopy} package.}

\references{
\itemize{
\item This data is copyright of the Soil Security Laboratory, The University of Sydney. If any part of this data is to be used in any publication or report, please provide a citation:
Soil Security Laboratory, 2016. R package for soil spectroscopy. The University of Sydney, Sydney, Australia.
}
}

\examples{

library(spectroscopy)

data(soilSpecDat)

#Histogram of all the clay data
hist(soilSpecDat$clay)

#Plot the vis-NIR spectra of the first row
plot(c(350:2500), soilSpecDat[1,5:ncol(soilSpecDat)], type = "l", col="red", 
xlab="nm", ylab="reflectance")


}
\keyword{datasets}
