\name{CAAGlobe_DistanceBetweenPoints}
\alias{CAAGlobe_DistanceBetweenPoints}
\title{
CAAGlobe_DistanceBetweenPoints
}
\description{
CAAGlobe_DistanceBetweenPoints
}
\usage{
CAAGlobe_DistanceBetweenPoints(GeographicalLatitude1, GeographicalLongitude1, GeographicalLatitude2, GeographicalLongitude2)
}
\arguments{
  \item{GeographicalLatitude1}{
GeographicalLatitude1 The latitude of the position of the first point in degrees.
}
  \item{GeographicalLongitude1}{
GeographicalLongitude1 The latitude of the position of the first point in degrees.

}
  \item{GeographicalLatitude2}{
GeographicalLatitude2 The latitude of the position of the second point in degrees.
}
  \item{GeographicalLongitude2}{
GeographicalLongitude2 The latitude of the position of the first point in degrees.

}
}
\details{
}
\value{
The shortest distance between two known points on the surface of the Earth expressed in kilometres.
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
CAAGlobe_DistanceBetweenPoints(50, 0, 50, 60);
CAAGlobe_DistanceBetweenPoints(50, 0, 50, 1);
}
\keyword{ Distance }

