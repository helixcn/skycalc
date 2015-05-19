\name{CAAPrecession_PrecessEcliptic}
\alias{CAAPrecession_PrecessEcliptic}
\title{
CAAPrecession_PrecessEcliptic
}
\description{
CAAPrecession_PrecessEcliptic
}
\usage{
CAAPrecession_PrecessEcliptic(Lambda, Beta, JD0, JD)
}
\arguments{
  \item{Lambda}{
Lambda The ecliptical longitude in degrees.
}
  \item{Beta}{
Beta The ecliptical latitude of the star in degrees.
}
  \item{JD0}{
JD0 The date in Dynamical time corresponding to the initial epoch.
}
  \item{JD}{ 
JD The date in Dynamical time corresponding to the final epoch.
}
}
\details{
}
\value{
Returns the converted ecliptic coordinates in a CAA2DCoordinate class. The x value in the class corresponds to the ecliptical longitude in degrees and the y value corresponds to the ecliptical latitude in degrees.  
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
CAAPrecession_PrecessEcliptic(Lambda = 134, Beta = 12, JD0 = 2418882, JD =  2455382)
}
\keyword{ Precession }

