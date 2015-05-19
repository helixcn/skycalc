\name{CAAMoslemCalendar_MoslemToJulian}
\alias{CAAMoslemCalendar_MoslemToJulian}
\title{
CAAMoslemCalendar_MoslemToJulian
}
\description{
CAAMoslemCalendar_MoslemToJulian
}
\usage{
CAAMoslemCalendar_MoslemToJulian(Year, Month, Day)
}
\arguments{
  \item{Year}{
Year The year in the Moslem Calendar to convert.
}
  \item{Month}{
Month The month of the year in the Moslem Calendar (1 for Muharram to 12 for Dhu l-Hijja).
}
  \item{Day}{
Day The day of the month in the Moslem Calendar.
}
}
\details{
}
\value{
A class containing

Year The year in the Julian Calendar.

Month The month of the year in the Julian Calendar  (1 for January to 12 for December).

Day The day of the month in the Julian Calendar.
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
CAAMoslemCalendar_MoslemToJulian(Year = 1431, Month = 2, Day = 12)
}
\keyword{ Calendar }

