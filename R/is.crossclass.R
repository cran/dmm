is.crossclass <-
function(x)
#  is.crossclass()  -  test if x is cross-class specific
{
  cclass <- !is.withinclass(x)
  return(cclass)
}
