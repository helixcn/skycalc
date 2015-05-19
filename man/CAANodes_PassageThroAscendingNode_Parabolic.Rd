\name{CAANodes_PassageThroAscendingNode_Parabolic}
\alias{CAANodes_PassageThroAscendingNode_Parabolic}
\title{
CAANodes_PassageThroAscendingNode_Parabolic
}
\description{
CAANodes_PassageThroAscendingNode_Parabolic
}
\usage{
CAANodes_PassageThroAscendingNode_Parabolic(elements_w, elements_T, elements_q)
}
\arguments{
  \item{elements_w}{
w The argument of the perihelion in degrees.
}
  \item{elements_T}{
T The Julian date of the time of passage in perihelion.
}
  \item{elements_q}{
q The perihelion distance in astronomical units.
}
}
\details{
}
\value{
A class containing

t The date in Dynamical time to when the body moves thro the ascending node.

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
CAANodes_PassageThroAscendingNode_Parabolic(elements_w =22, elements_T=2323553, elements_q=0.13)
}
\keyword{ Parabolic }

