is.specific <-
function(x)
#  is.specific()  -  test if x is name of a class-specific component
{
   colons <- which(strsplit(x,"")[[1]] == ":")
  if(length(colons) == 0) {  # nonspecific case
    specific <- F
  }
  else {
    specific <- T
  }
  return(specific)
}
