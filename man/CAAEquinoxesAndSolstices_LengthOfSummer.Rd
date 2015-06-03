\name{CAAEquinoxesAndSolstices_LengthOfSummer}
\alias{CAAEquinoxesAndSolstices_LengthOfSummer}

\title{
Length Of Summer
}
\description{
Compute the length of the astronomical Summer season in days.
}
\usage{
CAAEquinoxesAndSolstices_LengthOfSummer(Year, bNorthernHemisphere = TRUE)
}

\arguments{
  \item{Year}{
The year to calculate for.
}
  \item{bNorthernHemisphere}{
logical, true to indicate the observer is in the northern hemisphere, while false means the southern hemisphere.
}
}
\details{
The length of the astronomical Summer season in days.
}
\value{
The length of the astronomical Summer season in days.
}
\references{
    \url{http://www.naughter.com/aa.html}
}
\author{
 C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{

}


\seealso{

}
\examples{
    CAAEquinoxesAndSolstices_LengthOfSummer(Year = 2012, bNorthernHemisphere = TRUE)
}

\keyword{ season }

