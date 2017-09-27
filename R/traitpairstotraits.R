traitpairstotraits <-
function(tp)
# traitpairstotraits() - extract vector of traitnames from vector of traitpairs
{
  t <- rep(" ",length(tp))
  for(i in 1:length(tp)) {
    colon <- which(strsplit(tp[i],"")[[1]] == ":")
    t[i] <- substr(tp[i],1,colon-1)
  }
  tu <- unique(t)
  return(tu)
}
