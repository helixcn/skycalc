\name{CAAParallactic_ParallacticAngle}
\alias{CAAParallactic_ParallacticAngle}
\title{
CAAParallactic_ParallacticAngle
}
\description{
CAAParallactic_ParallacticAngle
}
\usage{
CAAParallactic_ParallacticAngle(HourAngle, Latitude, delta)
}
\arguments{
  \item{HourAngle}{
HourAngle The hour angle.
}
  \item{Latitude}{
Latitude The latitude of the position in degrees.
}
  \item{delta}{
Delta The declination in degrees.
}
}
\details{
}
\value{
Returns the parallactic angle (the ZCN angle) and is generally designated by q.
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
CAAParallactic_ParallacticAngle(HourAngle = 11.3, Latitude = 122, delta = 15)
}
\keyword{ Parallactic }

