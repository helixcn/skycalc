\name{CAAGalileanMoons_Calculate}
\alias{CAAGalileanMoons_Calculate}
\title{
CAAGalileanMoons_Calculate
}
\description{
CAAGalileanMoons_Calculate
}
\usage{
CAAGalileanMoons_Calculate(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 

}
}
\details{
}
\value{
A class which itself contains a class for each moon. This contained class itself contains

MeanLongitude The mean longitude of the moon in degrees.

TrueLongitude The true longitude of the moon in degrees.

TropicalLongitude The tropical longitude of the moon in degrees.

EquatorialLatitude The latitude in degrees of the moon with respect to Jupiter's equatorial plane.

r The radius vector of the moon in equatorial radii of Jupiter.

TrueRectangularCoordinates The true 3D rectangular coordinates of the moon.

ApparentRectangularCoordinates The apparent 3D rectangular coordinates of the moon.

bInTransit A Boolean which if true means that the moon is in front of Jupiter as viewed from the Earth otherwise false.

bInOccultation A Boolean which if true means that the moon is behind Jupiter as viewed from the Earth otherwise false.

bInEclipse A Boolean which if true means that the moon is behind Jupiter as viewed from the Sun otherwise false.

bInShadowTransit A Boolean which if true means that the moon is in front of Jupiter as viewed from the Earth otherwise false.
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
CAAGalileanMoons_Calculate(2448972.50068);
}
\keyword{ Jupiter }
\keyword{ Satellite }
