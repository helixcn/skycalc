\name{CAAPlanetPerihelionAphelion_MarsPerihelion}
\alias{CAAPlanetPerihelionAphelion_MarsPerihelion}
\title{
CAAPlanetPerihelionAphelion_MarsPerihelion
}
\description{
CAAPlanetPerihelionAphelion_MarsPerihelion
}
\usage{
CAAPlanetPerihelionAphelion_MarsPerihelion(k)
}
\arguments{
  \item{k}{
k The K value to calculate the perihelion for.
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
 options(digits = 15)
CAAPlanetPerihelionAphelion_MarsPerihelion(5)
}
\keyword{ Mars }
\keyword{ Perihelion }
