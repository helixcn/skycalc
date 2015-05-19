\name{CAAFK5_ConvertVSOPToFK5AnyEquinox}
\alias{CAAFK5_ConvertVSOPToFK5AnyEquinox}
\title{
CAAFK5_ConvertVSOPToFK5AnyEquinox
}
\description{
CAAFK5_ConvertVSOPToFK5AnyEquinox
}
\usage{
CAAFK5_ConvertVSOPToFK5AnyEquinox(value_X, value_Y, value_Z, JDEquinox)
}
\arguments{
  \item{value_X}{
value X The geometric rectangular ecliptical coordinates of the object (e.g. the Sun) to convert from the dynamical reference frame (VSOP) to the equatorial FK5 reference frame of JDEquinox.
}
  \item{value_Y}{
value Y The geometric rectangular ecliptical coordinates of the object (e.g. the Sun) to convert from the dynamical reference frame (VSOP) to the equatorial FK5 reference frame of JDEquinox.
}
  \item{value_Z}{
value Z The geometric rectangular ecliptical coordinates of the object (e.g. the Sun) to convert from the dynamical reference frame (VSOP) to the equatorial FK5 reference frame of JDEquinox.
}
  \item{JDEquinox}{
The Julian day for which equatorial coordinates should be calculated for.
}
}
\details{
}
\value{
A class containing the converted equatorial FK5 reference frame coordinates.
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
CAAFK5_ConvertVSOPToFK5AnyEquinox(value_X = 0.12, value_Y = 0.777, value_Z = 0.23, JDEquinox = 2456923)
}
\keyword{ FK5 }

