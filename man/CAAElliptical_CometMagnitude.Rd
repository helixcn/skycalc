\name{CAAElliptical_CometMagnitude}
\alias{CAAElliptical_CometMagnitude}
\title{
CAAElliptical_CometMagnitude
}
\description{
CAAElliptical_CometMagnitude
}
\usage{
CAAElliptical_CometMagnitude(g, delta, k, r)
}
\arguments{
  \item{g}{
g The absolute magnitude of the comet.
}
  \item{delta}{
delta Distance of the comet to the Earth in astronomical units.
}
  \item{k}{
k A constant which differs from one comet to another.
}
  \item{r}{
r Distance of the comet from the Sun in astronomical units
}
}
\details{
}
\value{
The magnitude of the comet.
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
Mag2 = CAAElliptical_CometMagnitude(5.5, 0.378, 10, 0.658);
Mag3 = CAAElliptical_CometMagnitude(5.5, 1.1017, 10, 1.5228);
}
\keyword{ Elliptical }

