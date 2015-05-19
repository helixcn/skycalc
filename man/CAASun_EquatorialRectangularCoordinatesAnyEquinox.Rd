\name{CAASun_EquatorialRectangularCoordinatesAnyEquinox}
\alias{CAASun_EquatorialRectangularCoordinatesAnyEquinox}
\title{
CAASun_EquatorialRectangularCoordinatesAnyEquinox
}
\description{
CAASun_EquatorialRectangularCoordinatesAnyEquinox
}
\usage{
CAASun_EquatorialRectangularCoordinatesAnyEquinox(JD, JDEquinox)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
  \item{JDEquinox}{
The Julian day for which equatorial coordinates should be calculated for.
}
}
\details{
}
\value{
A class containing the equatorial 3D rectangular coordinates in astronomical units referred to the specified equinox defined in the FK5 theory.
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
CAASun_EquatorialRectangularCoordinatesAnyEquinox(JD =  2456597.5, JDEquinox = 2000.0)
}
\keyword{ Sun }

