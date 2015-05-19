\name{CAADate_DaysInMonth}
\alias{CAADate_DaysInMonth}
\title{
CAADate_DaysInMonth
}
\description{
CAADate_DaysInMonth
}
\usage{
CAADate_DaysInMonth(Month, bLeap)
}
\arguments{
  \item{Month}{
Month: 1 to 12
}
  \item{bLeap}{
Logical, indicate whether the month is Leap or not. 
}
}
\details{
}
\value{
Returns the total number of days in the month (28 - 31) which this date represents. The static version of the function can be used if you do not want to construct a CAADate instance to do this test.
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
CAADate_DaysInMonth(Month = 2, bLeap = FALSE)
CAADate_DaysInMonth(Month = 2, bLeap = TRUE)
CAADate_DaysInMonth(Month = 3, bLeap = FALSE)
CAADate_DaysInMonth(Month = 3, bLeap = TRUE)
}
\keyword{ Date }

