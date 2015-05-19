\name{CAANutation_NutationInDeclination}
\alias{CAANutation_NutationInDeclination}
\title{
CAANutation_NutationInDeclination
}
\description{
CAANutation_NutationInDeclination
}
\usage{
CAANutation_NutationInDeclination(Alpha, Obliquity, NutationInLongitude, NutationInObliquity)
}
\arguments{
  \item{Alpha}{
Alpha The right ascension of the position in hour angles.
}
  \item{Obliquity}{
Obliquity The obliquity of the Ecliptic in degrees.
}
  \item{NutationInLongitude}{
NutationInLongitude The nutation in longitude in arc seconds of a degree.
}
  \item{NutationInObliquity}{
NutationInObliquity The nutation in obliquity in arc seconds of a degree.
}
}
\details{
}
\value{
The nutation in declination in arc seconds of a degree.
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
CAANutation_NutationInDeclination(Alpha = 12.33, 
Obliquity = 23.4444, NutationInLongitude = 0.124, 
NutationInObliquity = 0.22)
}
\keyword{ Nutation }

