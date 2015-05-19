\name{CAAJewishCalendar_IsLeap}
\alias{CAAJewishCalendar_IsLeap}
\title{
CAAJewishCalendar_IsLeap
}
\description{
CAAJewishCalendar_IsLeap
}
\usage{
CAAJewishCalendar_IsLeap(Year)
}
\arguments{
  \item{Year}{ Year The year in the Jewish calendar.
}
}
\details{
}
\value{
true if the specified year is leap otherwise false.
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
CAAJewishCalendar_IsLeap(Year = 1309)
CAAJewishCalendar_IsLeap(Year = 1305)
}
\keyword{ Calendar }

