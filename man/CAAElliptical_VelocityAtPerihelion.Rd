\name{CAAElliptical_VelocityAtPerihelion}
\alias{CAAElliptical_VelocityAtPerihelion}
\title{
CAAElliptical_VelocityAtPerihelion
}
\description{
CAAElliptical_VelocityAtPerihelion
}
\usage{
CAAElliptical_VelocityAtPerihelion(e, a)
}
\arguments{
  \item{e}{
e The eccentricity of the orbit.
}
  \item{a}{
a The semi major axis of the orbit in astronomical units.
}
}
\details{
}
\value{
The velocity of the object in kilometres per second at perihelion.
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
CAAElliptical_VelocityAtPerihelion(e = 0.21, a = 1.5)
}
\keyword{ Perihelion }

