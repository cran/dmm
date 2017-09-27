is.var <-
function(x)
#  is.var()   -   test if arg is rowname of a variance
# this may be rdundant - see is.cecov() and is.specific()
{
  colons <- which(strsplit(x,"")[[1]] == ":")
  if(length(colons) == 0) {  # nonspecific case
    isvar <- T
  }
  else {  # specific case
    class1 <- substr(x,colons[1]+1,colons[2]-1)
    class2 <- substr(x,colons[2]+1,colons[3]-1)

    if(class1 == class2) {
      isvar <- T
    }
    else {
      isvar <- F
    }
  }
  return(isvar)
}
