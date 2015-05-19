\name{CAASidereal_MeanGreenwichSiderealTime}
\alias{CAASidereal_MeanGreenwichSiderealTime}
\title{
CAASidereal_MeanGreenwichSiderealTime
}
\description{
CAASidereal_MeanGreenwichSiderealTime
}
\usage{
CAASidereal_MeanGreenwichSiderealTime(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{

}
\value{
The Mean Greenwich Sidereal Time, that is, the Greenwich hour angle of the mean vernal point (the intersection of the ecliptic of the date with the mean equator of the date), expressed in hours.
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
CAASidereal_MeanGreenwichSiderealTime(JD)
}
\keyword{ Time }

