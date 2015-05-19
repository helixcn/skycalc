\name{CAACoordinateTransformation_Equatorial2Horizontal}
\alias{CAACoordinateTransformation_Equatorial2Horizontal}
\title{
CAACoordinateTransformation_Equatorial2Horizontal
}
\description{
CAACoordinateTransformation_Equatorial2Horizontal
}
\usage{
CAACoordinateTransformation_Equatorial2Horizontal(LocalHourAngle, Delta, Latitude)
}
\arguments{
  \item{LocalHourAngle}{
LocalHourAngle The local hour angle, measured westwards from the South. 
}
  \item{Delta}{
Delta The declination in degrees.
}
  \item{Latitude}{
Latitude The standard latitude of the position in degrees.
}
}
\details{
}
\value{
Returns the converted horizontal coordinates in a CAA2DCoordinate class. The x value in the class corresponds to the azimuth in degrees and the y value corresponds to the altitude in degrees.
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
CAACoordinateTransformation_Equatorial2Horizontal(LocalHourAngle = 8.12, Delta = 11, Latitude = 15)
}
\keyword{ transformation }

