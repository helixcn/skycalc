\name{CAAPrecession_PrecessEquatorialFK4}
\alias{CAAPrecession_PrecessEquatorialFK4}
\title{
CAAPrecession_PrecessEquatorialFK4
}
\description{
CAAPrecession_PrecessEquatorialFK4
}
\usage{
CAAPrecession_PrecessEquatorialFK4(Alpha, Delta, JD0, JD)
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
CAAPrecession_PrecessEquatorialFK4(Alpha = 22, Delta = 41, JD0 = 285700, JD = 287900)
}
\keyword{ Precession }

