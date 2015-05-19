\name{CAASaturnMoons_Calculate}
\alias{CAASaturnMoons_Calculate}
\title{
CAASaturnMoons_Calculate
}
\description{
CAASaturnMoons_Calculate
}
\usage{
CAASaturnMoons_Calculate(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
A class which itself contains a class for each moon. This contained class itself contains

TrueRectangularCoordinates The true 3D rectangular coordinates of the moon.

ApparentRectangularCoordinates The apparent 3D rectangular coordinates of the moon.

bInTransit A Boolean which if true means that the moon is in front of Saturn as viewed from the Earth otherwise false.

bInOccultation A Boolean which if true means that the moon is behind Saturn as viewed from the Earth otherwise false.

bInEclipse A Boolean which if true means that the moon is behind Saturn as viewed from the Sun otherwise false.

bInShadowTransit A Boolean which if true means that the moon is in front of Saturn as viewed from the Earth otherwise false.
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
CAASaturnMoons_Calculate(JD)
}
\keyword{ Saturn }

