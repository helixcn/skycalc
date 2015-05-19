CAADate_DateToJD <-
function(Year, Month, Day, bGregorianCalendar){
res <- .Call("CAADate_DateToJD", Year, Month, Day, bGregorianCalendar)
return(res)
}
