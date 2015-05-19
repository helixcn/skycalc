\name{CAAElliptical_SemiMajorAxisFromPerihelionDistance}
\alias{CAAElliptical_SemiMajorAxisFromPerihelionDistance}
\title{
CAAElliptical_SemiMajorAxisFromPerihelionDistance
}
\description{
CAAElliptical_SemiMajorAxisFromPerihelionDistance
}
\usage{
CAAElliptical_SemiMajorAxisFromPerihelionDistance(q, e)
}
\arguments{
  \item{q}{
q The perihelion distance in astronomical units.
}
  \item{e}{
e The eccentricity of the orbit.
}
}
\details{
}
\value{
The semi major axis of the orbit in astronomical units.
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
CAAElliptical_SemiMajorAxisFromPerihelionDistance(q = 0.84, e = 0.12)
}
\keyword{ Perihelion }

