\name{CAAPhysicalJupiter_Calculate}
\alias{CAAPhysicalJupiter_Calculate}
\title{
CAAPhysicalJupiter_Calculate
}
\description{
CAAPhysicalJupiter_Calculate
}
\usage{
CAAPhysicalJupiter_Calculate(JD)
}
\arguments{
  \item{JD}{ JD The date in Dynamical time to calculate for.

}
}
\details{

}
\value{
A class containing

DE The planetocentric declination in degrees of the Earth.

DS The planetocentric declination in degrees of the Sun.

Geometricw1 The geometric longitude in degrees of the central meridian for System 1.

Geometricw2 The geometric longitude in degrees of the central meridian for System 2.

Apparentw1 The apparent longitude in degrees of the central meridian for System 1.

Apparentw2 The geometric longitude in degrees of the central meridian for System 2.

P The position angle of Jupiter's northern rotational pole in degrees.
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
CAAPhysicalJupiter_Calculate(JD  = 2456597.5) 
}
\keyword{ Jupiter }

