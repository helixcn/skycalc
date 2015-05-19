\name{CAACoordinateTransformation_Equatorial2Galactic}
\alias{CAACoordinateTransformation_Equatorial2Galactic}
\title{
CAACoordinateTransformation_Equatorial2Galactic
}
\description{
CAACoordinateTransformation_Equatorial2Galactic
}
\usage{
CAACoordinateTransformation_Equatorial2Galactic(Alpha, Delta)
}
\arguments{
  \item{Alpha}{
Alpha The right ascension expressed as an hour angle.
} 
  \item{Delta}{
Delta The declination in degrees.
}
}
\details{
Returns the converted galactic coordinates in a CAA2DCoordinate class. The x value in the class corresponds to the galactic longitude in degrees and the y value corresponds to the galactic latitude in degrees.
}
\value{
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
CAACoordinateTransformation_Equatorial2Galactic(Alpha = 12.3, Delta = 22)
}
\keyword{ transformation }

