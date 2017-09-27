match.vars <-
function(longrowname)
# match.vars()    - find 2 longnames for vars matching longrowname)
{
  colons <- which(strsplit(longrowname,"")[[1]] == ":")
  if(length(colons) == 0) {  # nonspecific case never a cros-class cov
    stop("match.vars() - should never get here:\n")
  }
  else { # specific case
    class1 <- substr(longrowname,colons[1]+1,colons[2]-1)
    class2 <- substr(longrowname,colons[2]+1,colons[3]-1)
    genericname <- substr(longrowname,colons[3]+1,nchar(longrowname))
    var1name <- paste(substr(longrowname,1,colons[1]-1),":",class1,":",class1,":",genericname,sep="")
    var2name <- paste(substr(longrowname,1,colons[1]-1),":",class2,":",class2,":",genericname,sep="")
    var1class <- paste(substr(longrowname,1,colons[1]-1),":",class1,":",class1,sep="")
    var2class <- paste(substr(longrowname,1,colons[1]-1),":",class2,":",class2,sep="")
  }
  outlist <- list(var1name=var1name,var2name=var2name,var1class=var1class,var2class=var2class,genericname=genericname)
  return(outlist)
}
