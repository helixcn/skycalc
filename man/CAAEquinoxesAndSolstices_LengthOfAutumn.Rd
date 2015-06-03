\name{CAAEquinoxesAndSolstices_LengthOfAutumn}
\alias{CAAEquinoxesAndSolstices_LengthOfAutumn}

\title{
The length of the astronomical Autumn season in days.
}
\description{
The length of the astronomical Autumn season in days.
}
\usage{
CAAEquinoxesAndSolstices_LengthOfAutumn(Year, bNorthernHemisphere = TRUE)
}

\arguments{
  \item{Year}{
Year The year to calculate the occurrence for. Note that this method refers to "Northward Equinox" instead of "Spring Equinox" to avoid a northern hemisphere-specific bias.
}
  \item{bNorthernHemisphere}{
Logical, implying whether the calculation is for northern hemisphere.
}
}
\details{
to be added. 
}
\value{
The length of the astronomical Autumn season in days.
}
\references{
    \url{http://www.naughter.com/aa.html}
}
\author{
 C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{

}

\seealso{

}
\examples{
CAAEquinoxesAndSolstices_LengthOfAutumn(Year = 2006, bNorthernHemisphere = TRUE)
}

\keyword{ season }
\keyword{ autumn }
