\name{CAARefraction_RefractionFromApparent}
\alias{CAARefraction_RefractionFromApparent}
\title{
CAARefraction_RefractionFromApparent
}
\description{
CAARefraction_RefractionFromApparent
}
\usage{
CAARefraction_RefractionFromApparent(Altitude, Pressure, Temperature)
}
\arguments{
  \item{Altitude}{
Altitude The apparent altitude in degrees.
}
  \item{Pressure}{
Pressure The atmospheric pressure in millibars
}
  \item{Temperature}{
Temperature The air temperature in degrees Celsius.
}
}
\details{
}
\value{
the refraction in degrees.
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
CAARefraction_RefractionFromApparent(25, 990, 25)
}
\keyword{ Refraction }

