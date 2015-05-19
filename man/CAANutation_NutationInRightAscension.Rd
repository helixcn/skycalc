\name{CAANutation_NutationInRightAscension}
\alias{CAANutation_NutationInRightAscension}
\title{
CAANutation_NutationInRightAscension
}
\description{
CAANutation_NutationInRightAscension
}
\usage{
CAANutation_NutationInRightAscension(Alpha, Delta, Obliquity, NutationInLongitude, NutationInObliquity)
}
\arguments{
  \item{Alpha}{
Alpha The right ascension of the position in hour angles.
}
  \item{Delta}{
Delta The declination of the position in degrees.
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
The nutation in right ascension in arc seconds of a degree.
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
CAANutation_NutationInRightAscension(Alpha = 10.93, Delta = 11.34, Obliquity = 23.54343, NutationInLongitude = 0.231, NutationInObliquity = 0.004)
}
\keyword{ Nutation }

