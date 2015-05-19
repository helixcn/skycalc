\name{CAAPlanetPerihelionAphelion_EarthPerihelion}
\alias{CAAPlanetPerihelionAphelion_EarthPerihelion}
\title{
CAAPlanetPerihelionAphelion_EarthPerihelion
}
\description{
CAAPlanetPerihelionAphelion_EarthPerihelion
}
\usage{
CAAPlanetPerihelionAphelion_EarthPerihelion(k, bBarycentric)
}
\arguments{
  \item{k}{
k The K value to calculate the perihelion for.
}
  \item{bBarycentric}{
bBarycentric  If true, the calculation is for the barycenter of the Earth - Moon system, false implies the Earth itself.
}
}
\details{
}
\value{
Returns the date in Dynamical time when the specified perihelion occurs.
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
CAAPlanetPerihelionAphelion_EarthPerihelion(k = 10, bBarycentric = TRUE)
CAAPlanetPerihelionAphelion_EarthPerihelion(k = 10, bBarycentric = FALSE)

}
\keyword{ Earth }
\keyword{ Perihelion }
