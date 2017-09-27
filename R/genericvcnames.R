genericvcnames <-
function(x)
#  genericvcnames()  - find the variance component part of a specific vc name
#                x is a vector of specific vc names
{
  gvcname <- rep(0,length(x))
  for( i in 1:length(x)) {
    colons <- which(strsplit(x[i],"")[[1]] == ":")
    if(length(colons) == 0) {
      gvcname[i] <- x[i]
    }
    else {
      lastcolon <- colons[length(colons)]
      gvcname[i] <- substr(x[i],lastcolon+1,nchar(x[i]))
    }
  }
  return(gvcname)
}
