\name{CAAPrecession_EquatorialPMToEcliptic}
\alias{CAAPrecession_EquatorialPMToEcliptic}
\title{
CAAPrecession_EquatorialPMToEcliptic
}
\description{
CAAPrecession_EquatorialPMToEcliptic
}
\usage{
CAAPrecession_EquatorialPMToEcliptic(Alpha, Delta, Beta, PMAlpha, PMDelta, Epsilon)
}
\arguments{
  \item{Alpha}{
Alpha The right ascension expressed as an hour angle.
}
  \item{Delta}{
Delta The declination in degrees.
}
  \item{Beta}{
Beta The ecliptical latitude of the star in degrees.
}
  \item{PMAlpha}{
PMAlpha The proper motion in right ascension in arc seconds per year.
}
  \item{PMDelta}{
PMDelta The proper motion of the declination in arc seconds per year.
}
  \item{Epsilon}{
Epsilon The obliquity of the ecliptic in degrees.
}
}
\details{
}
\value{
Returns the converted ecliptic proper motions in a CAA2DCoordinate class. The x value in the class corresponds to the proper motion in ecliptical longitude in arc seconds per year and the y value corresponds to the proper motion in ecliptical latitude in arc seconds per year.
}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated.
}
\author{ 
 C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
}
\seealso{
}
\examples{
CAAPrecession_EquatorialPMToEcliptic(Alpha = 12.54, Delta = 10, 
    Beta = 67, PMAlpha = 0.01, PMDelta = 0.21, Epsilon = 0.1)
}
\keyword{ Precession }

