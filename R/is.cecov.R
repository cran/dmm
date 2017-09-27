is.cecov <-
function(x,allcovs)
#  is.cecov()  -  test if x is name of a cross-effect covariance
#              -  if is a spsecific longname extracts genericname
{
  colons <- which(strsplit(x,"")[[1]] == ":")
  if(length(colons) == 0) {  # nonspecific case
    genericname <- x
  }
  else {
    genericname <- substr(x,colons[3]+1,nchar(x))
  }
  if(any(!is.na(match(genericname,allcovs)))) {
    cecov <- T
  }
  else {
    cecov <- F
  }
  return(cecov)
}
