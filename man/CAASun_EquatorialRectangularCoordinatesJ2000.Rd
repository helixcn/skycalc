\name{CAASun_EquatorialRectangularCoordinatesJ2000}
\alias{CAASun_EquatorialRectangularCoordinatesJ2000}
\title{
CAASun_EquatorialRectangularCoordinatesJ2000
}
\description{
CAASun_EquatorialRectangularCoordinatesJ2000
}
\usage{
CAASun_EquatorialRectangularCoordinatesJ2000(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
A class containing the equatorial 3D rectangular coordinates in astronomical units referred to the J2000 equinox defined in the FK5 theory.
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
CAASun_EquatorialRectangularCoordinatesJ2000(JD)
}
\keyword{ Sun }

