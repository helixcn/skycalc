\name{CAADynamicalTime_CumulativeLeapSeconds}
\alias{CAADynamicalTime_CumulativeLeapSeconds}

\title{
Compute Cumulative Leap Seconds
}
\description{
Compute Cumulative Leap Seconds
}
\usage{
CAADynamicalTime_CumulativeLeapSeconds(JD)
}

\arguments{
  \item{JD}{
Julian Date
}
}
\details{
See AAplus source code 1.60
}
\value{
the cumulative leap seconds
}
\references{
    \url{http://www.naughter.com/aa.html}
}
\author{
 C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
 This is a new function for AAplus v1.55 
}


\seealso{
    \code{\link{CAADynamicalTime_DeltaT}}
}

\examples{

CAADynamicalTime_CumulativeLeapSeconds(JD = 2456655.5)

}

\keyword{ time }

