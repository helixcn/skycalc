\name{CAARiseTransitSet_Calculate}
\alias{CAARiseTransitSet_Calculate}
\title{
CAARiseTransitSet_Calculate
}
\description{
CAARiseTransitSet_Calculate
}
\usage{
CAARiseTransitSet_Calculate(JD, Alpha1, Delta1, Alpha2, Delta2, Alpha3, Delta3, Longitude, Latitude, h0)
}
\arguments{
  \item{JD}{
  JD The Julian Day corresponding to that midnight Dynamical Time for the date when you want to perform the calculation.  
}
  \item{Alpha1}{ 
  Alpha1 The right ascension in hours of the object at time JD - 1 day
}
  \item{Delta1}{
Delta1 The declination in degrees of the object at time JD - 1 day
}
  \item{Alpha2}{
Alpha2 The right ascension in hours of the object at time JD
}
  \item{Delta2}{
Delta2 The declination in degrees of the object at time JD
}
  \item{Alpha3}{
Alpha3 The right ascension in hours of the object at time JD + 1 day
}
  \item{Delta3}{
Delta3 The declination in degrees of the object at time JD + 1 day
}
  \item{Longitude}{
Longitude The geographic longitude of the observer in degrees.
}
  \item{Latitude}{
Latitude The geographic latitude of the observer in degrees.
}
  \item{h0}{
h0 The "standard" altitude in degrees i.e. the geometric altitude of the centre of the body at the time of the apparent rising or setting. For stars and planets, you would normally use -0.5667, for the Sun you would use -0.8333 and for the moon you would use 0.7275 * PI - 0.5666 where PI is the Moon's horizontal parallax in degrees (If no great accuracy is required, the mean value of h0 = 0.125 can be used).
}
}
\details{
}
\value{

}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated.
}
\author{  C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
Ensure the RA values are corrected for interpolation. Due to important Remark 2 by Meeus on Interopolation of RA values
}
\seealso{
}
\examples{
CAARiseTransitSet_Calculate(JD = 2456597.5, Alpha1 = 123.1, Delta1 = 22.1 , 
Alpha2 = 124.2, Delta2 = 22.2, Alpha3 = 124.9, Delta3 = 22.3, 
Longitude = 117.20, Latitude = 46.45, h0 = 0.5667)
}
\keyword{ Rise }
\keyword{ Set }
\keyword{ skyevent }