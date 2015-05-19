\name{CAAEaster_Calculate}
\alias{CAAEaster_Calculate}
\title{
CAAEaster_Calculate
}
\description{
CAAEaster_Calculate
}
\usage{
CAAEaster_Calculate(nYear, GregorianCalendar)
}
\arguments{
  \item{nYear}{
nYear The year to perform the calculation for.
}
  \item{GregorianCalendar}{
bGregorianCalendar true if the calculation is to be performed for the Gregorian calendar, false implies the Julian Calendar
}
}
\details{
}
\value{
A class containing

Month The month on which Easter Sunday occurs.

Day The day of the month on which Easter Sunday occurs.
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
CAAEaster_Calculate(nYear = 1998, GregorianCalendar = TRUE)
CAAEaster_Calculate(nYear = 2014, GregorianCalendar = TRUE)
}
\keyword{ Easter }

