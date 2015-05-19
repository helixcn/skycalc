\name{CAASun_EquatorialRectangularCoordinatesMeanEquinox}
\alias{CAASun_EquatorialRectangularCoordinatesMeanEquinox}
\title{
CAASun_EquatorialRectangularCoordinatesMeanEquinox
}
\description{
CAASun_EquatorialRectangularCoordinatesMeanEquinox
}
\usage{
CAASun_EquatorialRectangularCoordinatesMeanEquinox(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
A class containing the equatorial 3D rectangular coordinates in astronomical units referred to the mean dynamical ecliptic and equinox of the date defined in the FK5 theory.
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
CAASun_EquatorialRectangularCoordinatesMeanEquinox(JD)
}
\keyword{ Sun }

