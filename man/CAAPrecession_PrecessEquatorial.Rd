\name{CAAPrecession_PrecessEquatorial}
\alias{CAAPrecession_PrecessEquatorial}
\title{
CAAPrecession_PrecessEquatorial
}
\description{
CAAPrecession_PrecessEquatorial
}
\usage{
CAAPrecession_PrecessEquatorial(Alpha, Delta, JD0, JD)
}
\arguments{
  \item{Alpha}{
Alpha The right ascension in hours of the object at time JD.
}
  \item{Delta}{
Delta The declination in degrees of the object at time JD.
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
CAAPrecession_PrecessEquatorial(Alpha = 35, Delta = 10, JD0 = 243211, JD = 264211)
}
\keyword{ Precession }

