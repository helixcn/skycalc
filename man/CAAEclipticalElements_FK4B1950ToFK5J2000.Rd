\name{CAAEclipticalElements_FK4B1950ToFK5J2000}
\alias{CAAEclipticalElements_FK4B1950ToFK5J2000}
\title{
CAAEclipticalElements_FK4B1950ToFK5J2000
}
\description{
CAAEclipticalElements_FK4B1950ToFK5J2000
}
\usage{
CAAEclipticalElements_FK4B1950ToFK5J2000(i0, w0, omega0)
}
\arguments{
  \item{i0}{
i0 The inclination in degrees to reduce .
}
  \item{w0}{
w0 The argument of perihelion in degrees to reduce.
}
  \item{omega0}{
omega0 The longitude of the ascending node in degrees to reduce.
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
CAAEclipticalElements_FK4B1950ToFK5J2000(i0 = 5.61234, w0 = 122.045, omega0 = 32.244)
}
\keyword{ FK4B1950}
\keyword{ FK5J2000 }
