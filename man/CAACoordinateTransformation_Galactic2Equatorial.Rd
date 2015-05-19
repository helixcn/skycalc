\name{CAACoordinateTransformation_Galactic2Equatorial}
\alias{CAACoordinateTransformation_Galactic2Equatorial}
\title{
CAACoordinateTransformation_Galactic2Equatorial
}
\description{
CAACoordinateTransformation_Galactic2Equatorial
}
\usage{
CAACoordinateTransformation_Galactic2Equatorial(l, b)
}
\arguments{
  \item{l}{
l  The galactic longitude expressed in degrees.
}
  \item{b}{
b The galactic latitude expressed in degrees.
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
CAACoordinateTransformation_Galactic2Equatorial(l = 35, b = 41)
}
\keyword{ transformation }

