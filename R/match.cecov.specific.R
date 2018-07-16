match.cecov.specific <-
function(longrowname)
# match.cecov.specific()    - find 2 longnames for vars matching a crosseffect
#                               longrowname
#                           - case 4  specific cross-class cross-effect cov
{
  colons <- which(strsplit(longrowname,"")[[1]] == ":")
  if(length(colons) == 0) {  # nonspecific case never a cross-class cov
    stop("match.cecov.specific() - should never get here:\n")
  }
  else { # specific case
    class1 <- substr(longrowname,colons[1]+1,colons[2]-1)
    class2 <- substr(longrowname,colons[2]+1,colons[3]-1)
    genericname <- substr(longrowname,colons[3]+1,nchar(longrowname))
    lbracket <- which(strsplit(genericname,"")[[1]] == "(")
    rbracket <- which(strsplit(genericname,"")[[1]] == ")")
    comma <- which(strsplit(genericname,"")[[1]] == ",")
    eff1 <- substr(genericname,lbracket[1] + 1, comma[1] - 1)
    eff2 <- substr(genericname,comma[1] + 1, rbracket[1]-1)
    genericname1 <- paste("Var",substr(genericname,4,4),"(",eff1,")",sep="")
    genericname2 <- paste("Var",substr(genericname,4,4),"(",eff2,")",sep="")
    var1name <- paste(substr(longrowname,1,colons[1]-1),":",class1,":",class1,":",genericname1,sep="")
    var2name <- paste(substr(longrowname,1,colons[1]-1),":",class2,":",class2,":",genericname2,sep="")
    var1class <- paste(substr(longrowname,1,colons[1]-1),":",class1,":",class1,sep="")
    var2class <- paste(substr(longrowname,1,colons[1]-1),":",class2,":",class2,sep="")
  }
  outlist <- list(var1name=var1name,var2name=var2name,var1class=var1class,var2class=var2class,genericname=genericname,genericname1=genericname1,genericname2=genericname2)
  return(outlist)
}
