\name{CAAEclipticalElements_Calculate}
\alias{CAAEclipticalElements_Calculate}
\title{
CAAEclipticalElements_Calculate
}
\description{
CAAEclipticalElements_Calculate
}
\usage{
CAAEclipticalElements_Calculate(i0, w0, omega0, JD0, JD)
}
\arguments{
  \item{i0}{
i0 The inclination in degrees to reduce.
}
  \item{w0}{
w0 The argument of perihelion in degrees to reduce.
}
  \item{omega0}{ 
omega0 The longitude of the ascending node in degrees to reduce.
}
  \item{JD0}{
JD0 The initial epoch in Dynamical time to calculate for.
}
  \item{JD}{ The date in Dynamical time to calculate for.
JD The epoch in Dynamical time to reduce the elements to.
}
}
\details{
}
\value{
A struct containing the following values:
i The reduced inclination in degrees.
w The reduced argument of perihelion in degrees.
omega The reduced longitude of the ascending node in degrees
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
CAAEclipticalElements_Calculate(i0 = 6.322, w0 = 34.23, 
omega0 = 45.11, JD0 = 2454555.5, JD= 2456623.5)
}
\keyword{ Ecliptical }
\keyword{ orbit }
