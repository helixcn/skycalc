CAAElliptical_Calculate <-
function(JD, elements_a, elements_e, elements_i, elements_w, elements_omega, elements_JDEquinox, elements_T){
.Call("CAAElliptical_Calculate", JD, elements_a, elements_e, elements_i, elements_w, elements_omega, elements_JDEquinox, elements_T)
}
