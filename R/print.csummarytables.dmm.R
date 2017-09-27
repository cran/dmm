print.csummarytables.dmm <-
function(x,...)
# print.csummarytables.dmm()  -  format a csummarytables.dmm object for printing
{
  for(i in 1: length(x$ctables)) {
    print(x$ctables[[i]],digits=x$digits)
    cat("\n")
  }
}
