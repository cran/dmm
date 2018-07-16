match.cecov.nonspecific <-
function(longrowname)
# match.cecov.nonspecific()   - find 2 longnames for vars matching a crosseffect#                               longrowname
#                             - case 3 - nonspecific cross-effect cov
{
  genericname <- substr(longrowname, 1 ,nchar(longrowname))
    lbracket <- which(strsplit(genericname,"")[[1]] == "(")
    rbracket <- which(strsplit(genericname,"")[[1]] == ")")
    comma <- which(strsplit(genericname,"")[[1]] == ",")
    eff1 <- substr(genericname,lbracket[1] + 1, comma[1] - 1)
    eff2 <- substr(genericname,comma[1] + 1, rbracket[1]-1)
    genericname1 <- paste("Var",substr(genericname,4,lbracket[1]-1),"(",eff1,")",sep="")
    genericname2 <- paste("Var",substr(genericname,4,lbracket[1]-1),"(",eff2,")",sep="")

    var1name <- genericname1
    var2name <- genericname2
    var1class <- NULL
    var2class <- NULL

  outlist <- list(var1name=var1name,var2name=var2name,var1class=var1class,var2class=var2class,genericname=genericname,genericname1=genericname1,genericname2=genericname2)
  return(outlist)
}
