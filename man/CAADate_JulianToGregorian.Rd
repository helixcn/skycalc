\name{CAADate_JulianToGregorian}
\alias{CAADate_JulianToGregorian}
\title{
CAADate_JulianToGregorian
}
\description{
CAADate_JulianToGregorian
}
\usage{
CAADate_JulianToGregorian(Year, Month, Day)
}
\arguments{
  \item{Year}{ Year The year. (Years are counted astronomically i.e. 1 BC = Year 0)
Year The year in the Julian Calendar to convert. (Years are counted astronomically i.e. 1 BC = Year 0)
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
Year The year in the Gregorian Calendar. (Years are counted astronomically i.e. 1 BC = Year 0)
Month The month of the year in the Gregorian Calendar (1 for January to 12 for December).
Day The day of the month in the Gregorian Calendar.
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
CAADate_JulianToGregorian(Year = 1999, Month = 9 , Day = 1.53312)
}
\keyword{ Date }

