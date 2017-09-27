is.withinclass <-
function(x)
# is.withinclass()  -  test if x is within-class ( assumes is specific)
{
  colons <- which(strsplit(x,"")[[1]] == ":") 
  if(length(colons) == 0 ) {  # nonspecific case
    wclass <- T
  }
  else {  # specific case
    class1 <- substr(x,colons[1]+1,colons[2]-1)
    class2 <- substr(x,colons[2]+1,colons[3]-1)
    if ( class1 == class2) {
      wclass <- T
    } 
    else {
      wclass <- F
    }
  }
  return(wclass)
}
