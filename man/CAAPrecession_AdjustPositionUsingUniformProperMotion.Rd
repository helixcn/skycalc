\name{CAAPrecession_AdjustPositionUsingUniformProperMotion}
\alias{CAAPrecession_AdjustPositionUsingUniformProperMotion}
\title{
CAAPrecession_AdjustPositionUsingUniformProperMotion
}
\description{
CAAPrecession_AdjustPositionUsingUniformProperMotion
}
\usage{
CAAPrecession_AdjustPositionUsingUniformProperMotion(t, Alpha, Delta, PMAlpha, PMDelta)
}
\arguments{
  \item{t}{
t The number of years from the starting epoch, negative in the past and positive in the future.
}
  \item{Alpha}{
Alpha The right ascension expressed as an hour angle.
}
  \item{Delta}{
Delta The declination in degrees.
}
  \item{PMAlpha}{
PMAlpha The proper motion in right ascension in arc seconds per year.
}
  \item{PMDelta}{
PMDelta The proper motion in declination in arc seconds per year.
}
}
\details{
}
\value{
Returns the converted equatorial coordinates in a CAA2DCoordinate class. The x value in the class corresponds to the right ascension expressed as an hour angle and the y value corresponds to the declination in degrees.
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
CAAPrecession_AdjustPositionUsingUniformProperMotion(t = 55, Alpha = 11.20, Delta = 45, PMAlpha = 1.22, PMDelta = 1.01)
}
\keyword{ Precession }

