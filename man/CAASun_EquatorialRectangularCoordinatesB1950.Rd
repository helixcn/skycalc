\name{CAASun_EquatorialRectangularCoordinatesB1950}
\alias{CAASun_EquatorialRectangularCoordinatesB1950}
\title{
CAASun_EquatorialRectangularCoordinatesB1950
}
\description{
CAASun_EquatorialRectangularCoordinatesB1950
}
\usage{
CAASun_EquatorialRectangularCoordinatesB1950(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{ 
A class containing the equatorial 3D rectangular coordinates in astronomical units referred to the B1950 equinox defined in the FK5 theory.
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
JD = 2456597.5
CAASun_EquatorialRectangularCoordinatesB1950(JD)
}
\keyword{ Sun }

