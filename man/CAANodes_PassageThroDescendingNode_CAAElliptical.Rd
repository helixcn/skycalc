\name{CAANodes_PassageThroDescendingNode_CAAElliptical}
\alias{CAANodes_PassageThroDescendingNode_CAAElliptical}
\title{
CAANodes_PassageThroDescendingNode_CAAElliptical
}
\description{
CAANodes_PassageThroDescendingNode_CAAElliptical
}
\usage{
CAANodes_PassageThroDescendingNode_CAAElliptical(elements_w, elements_e, elements_a, elements_T)
}
\arguments{
  \item{elements_w}{
w The argument of the perihelion in degrees.
}
  \item{elements_e}{
e The eccentricity of the orbit.
}
  \item{elements_a}{
a The semi major axis in astronomical units.
}
  \item{elements_T}{
T The Julian date of the time of passage in perihelion.
}
}
\details{
}
\value{
t The date in Dynamical time to when the body moves through the ascending node.

radius The radius vector in astronomical units.
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
CAANodes_PassageThroDescendingNode_CAAElliptical(elements_w = 113.1235, 
elements_e = 0.45, elements_a = 15.2, elements_T = 2457397.3434)
}
\keyword{ Elliptical }

