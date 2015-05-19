\name{CAAFK5_ConvertVSOPToFK5B1950}
\alias{CAAFK5_ConvertVSOPToFK5B1950}
\title{
CAAFK5_ConvertVSOPToFK5B1950
}
\description{
CAAFK5_ConvertVSOPToFK5B1950
}
\usage{
CAAFK5_ConvertVSOPToFK5B1950(value_X, value_Y, value_Z)
}
\arguments{
  \item{value_X}{
value X The geometric rectangular ecliptical coordinates of the object (e.g. the Sun) to convert from the dynamical reference frame (VSOP) to the equatorial FK5 B1950 reference frame.

}
  \item{value_Y}{
value Y The geometric rectangular ecliptical coordinates of the object (e.g. the Sun) to convert from the dynamical reference frame (VSOP) to the equatorial FK5 B1950 reference frame.

}
  \item{value_Z}{
value Z The geometric rectangular ecliptical coordinates of the object (e.g. the Sun) to convert from the dynamical reference frame (VSOP) to the equatorial FK5 B1950 reference frame.

}
}
\details{
}
\value{
A class containing the converted equatorial FK5 B1950 reference frame coordinates.
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
CAAFK5_ConvertVSOPToFK5B1950(value_X = 0.12, value_Y =0.34 , value_Z = 0.98)
}
\keyword{ FK5 }
\keyword{ transformation }
