\name{CAACoordinateTransformation_Ecliptic2Equatorial}
\alias{CAACoordinateTransformation_Ecliptic2Equatorial}
\title{
CAACoordinateTransformation_Ecliptic2Equatorial
}
\description{
CAACoordinateTransformation_Ecliptic2Equatorial
}
\usage{
CAACoordinateTransformation_Ecliptic2Equatorial(Lambda, Beta, Epsilon)
}
\arguments{
  \item{Lambda}{
Lambda The ecliptic longitude in degrees.
}
  \item{Beta}{
Beta The ecliptic latitude in degrees.
}
  \item{Epsilon}{
Epsilon the obliquity of the ecliptic in degrees.
}
}
\details{
}
\value{
Returns the converted equatorial coordinates in a CAA2DCoordinate class. The x value in the class corresponds to the equatorial longitude in hours and the y value corresponds to the equatorial latitude in degrees.
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
CAACoordinateTransformation_Ecliptic2Equatorial(Lambda = 34, Beta = 21, Epsilon = 0.0023)
}
\keyword{ transformation }

