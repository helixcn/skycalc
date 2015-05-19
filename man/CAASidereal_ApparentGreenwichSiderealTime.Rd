\name{CAASidereal_ApparentGreenwichSiderealTime}
\alias{CAASidereal_ApparentGreenwichSiderealTime}
\title{
CAASidereal_ApparentGreenwichSiderealTime
}
\description{
CAASidereal_ApparentGreenwichSiderealTime
}
\usage{
CAASidereal_ApparentGreenwichSiderealTime(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
The Apparent Greenwich Sidereal Time, that is, the Greenwich hour angle of the true vernal equinox, expressed in hours.
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
CAASidereal_ApparentGreenwichSiderealTime(JD)
}
\keyword{ Time }

