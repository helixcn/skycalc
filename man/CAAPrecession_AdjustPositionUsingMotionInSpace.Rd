\name{CAAPrecession_AdjustPositionUsingMotionInSpace}
\alias{CAAPrecession_AdjustPositionUsingMotionInSpace}
\title{
CAAPrecession_AdjustPositionUsingMotionInSpace
}
\description{
CAAPrecession_AdjustPositionUsingMotionInSpace
}
\usage{
CAAPrecession_AdjustPositionUsingMotionInSpace(r, DeltaR, t, Alpha, Delta, PMAlpha, PMDelta)
}
\arguments{
  \item{r}{
r The stars distance in parsecs.
}
  \item{DeltaR}{
deltar The radial velocity  in km/s.
}
  \item{t}{
t The number of years from the starting epoch, negative in the past and positive in the future.
}
  \item{Alpha}{
Alpha The right ascension expression as an hour angle.
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
CAAPrecession_AdjustPositionUsingMotionInSpace(r = 12, DeltaR = 250, t = 100,
    Alpha = 12.6, Delta = 22, PMAlpha = 1.2, PMDelta = 0.75)
}
\keyword{ Precession }

