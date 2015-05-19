\name{CAAMoslemCalendar_JulianToMoslem}
\alias{CAAMoslemCalendar_JulianToMoslem}
\title{
CAAMoslemCalendar_JulianToMoslem
}
\description{
CAAMoslemCalendar_JulianToMoslem
}
\usage{
CAAMoslemCalendar_JulianToMoslem(Year, Month, Day)
}
\arguments{
  \item{Year}{ Year The year. (Years are counted astronomically i.e. 1 BC = Year 0)
}
  \item{Month}{
Month The month of the year in the Julian Calendar (1 for January to 12 for December).
}
  \item{Day}{
Day The day of the month in the Julian Calendar.
}
}
\details{
}
\value{
A class containing

Year The year in the Moslem Calendar.

Month The month of the year in the Moslem Calendar (1 for Muharram to 12 for Dhu l-Hijja).

Day The day of the month in the Moslem Calendar.
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
CAAMoslemCalendar_JulianToMoslem(Year = 1999, Month = 1, Day = 1)
}
\keyword{ Calendar }

