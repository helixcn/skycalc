\name{CAAPlanetPerihelionAphelion_EarthAphelion}
\alias{CAAPlanetPerihelionAphelion_EarthAphelion}
\title{
CAAPlanetPerihelionAphelion_EarthAphelion
}
\description{
CAAPlanetPerihelionAphelion_EarthAphelion
}
\usage{
CAAPlanetPerihelionAphelion_EarthAphelion(k, bBarycentric)
}
\arguments{
  \item{k}{
k The K value to calculate the aphelion for.
}
  \item{bBarycentric}{
bBarycentric  If true, the calculation is for the barycenter of the Earth - Moon system, false implies the Earth itself
}
}
\details{
Returns the date in Dynamical time when the specified aphelion occurs.
}
\value{
}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated.
}
\author{ 
 C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
}
\seealso{
}
\examples{
CAAPlanetPerihelionAphelion_EarthAphelion(k = 10, bBarycentric = TRUE)
CAAPlanetPerihelionAphelion_EarthAphelion(k = 10, bBarycentric = FALSE)

}
\keyword{ Earth }
\keyword{ Aphelion }
