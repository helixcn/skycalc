CAAParabolic_Calculate <-
function(JD, elements_JDEquinox, elements_omega, elements_w, elements_i, elements_T, elements_q){
.Call("CAAParabolic_Calculate", JD, elements_JDEquinox, elements_omega, elements_w, elements_i, elements_T, elements_q)
}
