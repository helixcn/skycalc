\name{CAAParallactic_AngleBetweenNorthCelestialPoleAndNorthPoleOfEcliptic}
\alias{CAAParallactic_AngleBetweenNorthCelestialPoleAndNorthPoleOfEcliptic}
\title{
CAAParallactic_AngleBetweenNorthCelestialPoleAndNorthPoleOfEcliptic
}
\description{
CAAParallactic_AngleBetweenNorthCelestialPoleAndNorthPoleOfEcliptic
}
\usage{
CAAParallactic_AngleBetweenNorthCelestialPoleAndNorthPoleOfEcliptic(Lambda, Beta, ObliquityOfEcliptic)
}
\arguments{
  \item{Lambda}{
Lambda The ecliptical longitude in degrees.
}
  \item{Beta}{
Beta The ecliptical latitude of the star in degrees.
}
  \item{ObliquityOfEcliptic}{
ObliquityOfEcliptic The obliquity of the ecliptic in degrees.
}
}
\details{
}
\value{
Return the angle in degrees between the direction of the northern celestial pole and the direction of the north pole of the ecliptic, at the star.
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
CAAParallactic_AngleBetweenNorthCelestialPoleAndNorthPoleOfEcliptic(Lambda = 11, Beta = 45, ObliquityOfEcliptic = 33)
}
\keyword{ Parallactic }

