\name{CAACoordinateTransformation_Horizontal2Equatorial}
\alias{CAACoordinateTransformation_Horizontal2Equatorial}
\title{
CAACoordinateTransformation_Horizontal2Equatorial
}
\description{
CAACoordinateTransformation_Horizontal2Equatorial
}
\usage{
CAACoordinateTransformation_Horizontal2Equatorial(Azimuth, Altitude, Latitude)
}
\arguments{
  \item{Azimuth}{
A The azimuth in degrees.
}
  \item{Altitude}{
h The altitude in degrees
}
  \item{Latitude}{
Latitude The standard latitude of the position in degrees.
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
CAACoordinateTransformation_Horizontal2Equatorial(Azimuth = 110, Altitude = 35, Latitude = 41)
}
\keyword{ transformation }

