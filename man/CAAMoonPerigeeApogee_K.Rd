\name{CAAMoonPerigeeApogee_K}
\alias{CAAMoonPerigeeApogee_K}
\title{
CAAMoonPerigeeApogee_K
}
\description{
CAAMoonPerigeeApogee_K
}
\usage{
CAAMoonPerigeeApogee_K(Year)
}
\arguments{
  \item{Year}{ 
  Year The Year including decimals to calculate the K value for.
}
}
\details{
}
\value{
Returns the approximate value of K (required by the other methods of CAAMoonPerigeeApogee) for calculation of the specified phase.
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
CAAMoonPerigeeApogee_K(Year = 2015)
}
\keyword{ Moon }

