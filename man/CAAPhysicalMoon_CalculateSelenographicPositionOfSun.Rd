\name{CAAPhysicalMoon_CalculateSelenographicPositionOfSun}
\alias{CAAPhysicalMoon_CalculateSelenographicPositionOfSun}
\title{
CAAPhysicalMoon_CalculateSelenographicPositionOfSun
}
\description{
CAAPhysicalMoon_CalculateSelenographicPositionOfSun
}
\usage{
CAAPhysicalMoon_CalculateSelenographicPositionOfSun(JD)
}
\arguments{
  \item{JD}{ 
  The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
A class containing

l0 The longitude in degrees of the sub solar point.

b0 The latitude in degrees of the sub solar point.

c0 The selenographic colongitude of the sun.

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
CAAPhysicalMoon_CalculateSelenographicPositionOfSun(JD = 2455371.86)
}
\keyword{ Moon }
\keyword{ Sun }
