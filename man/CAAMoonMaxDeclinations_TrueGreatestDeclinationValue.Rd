\name{CAAMoonMaxDeclinations_TrueGreatestDeclinationValue}
\alias{CAAMoonMaxDeclinations_TrueGreatestDeclinationValue}
\title{
CAAMoonMaxDeclinations_TrueGreatestDeclinationValue
}
\description{
CAAMoonMaxDeclinations_TrueGreatestDeclinationValue
}
\usage{
CAAMoonMaxDeclinations_TrueGreatestDeclinationValue(k, bNortherly)
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
CAAMoonMaxDeclinations_TrueGreatestDeclinationValue(k = 45, bNortherly = TRUE)
}
\keyword{ Moon }

