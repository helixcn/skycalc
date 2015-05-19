\name{CAACoordinateTransformation_Equatorial2Ecliptic}
\alias{CAACoordinateTransformation_Equatorial2Ecliptic}
\title{
CAACoordinateTransformation_Equatorial2Ecliptic
}
\description{
CAACoordinateTransformation_Equatorial2Ecliptic
}
\usage{
CAACoordinateTransformation_Equatorial2Ecliptic(Alpha, Delta, Epsilon)
}
\arguments{
  \item{Alpha}{
Alpha The right ascension expressed as an hour angle. 
}
  \item{Delta}{
Delta The declination in degrees.
}
  \item{Epsilon}{
Epsilon The obliquity of the ecliptic in degrees.
}
}
\details{
}
\value{
Returns the converted ecliptic coordinates in a CAA2DCoordinate class. The x value in the class corresponds to the ecliptic longitude in degrees and the y value corresponds to the ecliptic latitude in degrees.
}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated. 
The transformation of coordinates from Equatorial to Ecliptic. This refers to algorithm 13.1 and 13.2 on page 93.
}
\author{  C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
}
\seealso{

}
\examples{
CAACoordinateTransformation_Equatorial2Ecliptic(Alpha = 12, Delta = 36, Epsilon = 23.33312)
}
\keyword{ transformation }

