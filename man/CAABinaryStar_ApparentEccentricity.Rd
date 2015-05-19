\name{CAABinaryStar_ApparentEccentricity}
\alias{CAABinaryStar_ApparentEccentricity}
\title{
CAABinaryStar_ApparentEccentricity
}
\description{
CAABinaryStar_ApparentEccentricity
}
\usage{
CAABinaryStar_ApparentEccentricity(e, i, w)
}
\arguments{
  \item{e}{
e Eccentricity of the true orbit.
}
  \item{i}{
i The inclination of the plane of the true orbit in degrees to the plane at right angles to the line of sight.
}
  \item{w}{
w The longitude of the periastron in degrees.
}
}
\details{
}
\value{
The apparent eccentricity of the binary star orbit 
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
CAABinaryStar_ApparentEccentricity(e = 0.70, i = 22, w = 123)
}
\keyword{ Binary Star }

