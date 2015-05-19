CAANearParabolic_Calculate <-
function(JD, elements_q, elements_i, elements_w, elements_omega, elements_JDEquinox, elements_T, elements_e){
.Call("CAANearParabolic_Calculate", JD, elements_q, elements_i, elements_w, elements_omega, elements_JDEquinox, elements_T, elements_e)
}
