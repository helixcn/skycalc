\name{CAAIlluminatedFraction_JupiterMagnitudeAA}
\alias{CAAIlluminatedFraction_JupiterMagnitudeAA}
\title{
CAAIlluminatedFraction_JupiterMagnitudeAA
}
\description{
CAAIlluminatedFraction_JupiterMagnitudeAA
}
\usage{
CAAIlluminatedFraction_JupiterMagnitudeAA(r, Delta, i)
}
\arguments{
  \item{r}{
r The planet's distance to the Sun in astronomical units.
}
  \item{Delta}{
Delta The planet's distance from the Earth in astronomical units.
}
  \item{i}{
i The planet's phase angle in degrees.
}
}
\details{
}
\value{
The magnitude of the planet.
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
CAAIlluminatedFraction_JupiterMagnitudeAA(r  = 12, Delta = 11.5,  i = 98)
}
\keyword{ Jupiter }

