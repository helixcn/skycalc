\name{CAAMoonMaxDeclinations_TrueGreatestDeclination}
\alias{CAAMoonMaxDeclinations_TrueGreatestDeclination}
\title{
CAAMoonMaxDeclinations_TrueGreatestDeclination
}
\description{
CAAMoonMaxDeclinations_TrueGreatestDeclination
}
\usage{
CAAMoonMaxDeclinations_TrueGreatestDeclination(k, bNortherly)
}
\arguments{
  \item{k}{
k The K value to calculate the max declination for.
}
  \item{bNortherly}{
bNortherly true if this is a calculation for a maximum northerly declination, false implies a calculation of the a maximum southerly declination.
}
}
\details{
}
\value{
Returns the date in Dynamical time when the true maximum declination occurs.
}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated.
}
\author{  C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
}
\seealso{
}
\examples{
CAAMoonMaxDeclinations_TrueGreatestDeclination(k = 45, bNortherly = TRUE)
}
\keyword{ Moon }

