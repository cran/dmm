make.csummarytables <-
function(object,traitset="all",componentset="all",bytrait=T,gls=F,digits=3)
# make.csummarytables()  -  make one set of component summary tables
{
    if(traitset[1] == "all"){
      traitpairs <- dimnames(object$variance.components)[[2]]
      traits <- traitpairstotraits(traitpairs)
    }
    else {
      traits <- traitset
      traitpairs <- permpaste(traits)
    }
    l <- length(traits)
    alltraitpairs <- dimnames(object$variance.components)[[2]]
 
    if(componentset[1] == "all") {
      components <- dimnames(object$variance.components)[[1]]
    }
    else {
      components <- componentset
    }

  if(bytrait) {

    ctables <- vector("list",l*l)  # one table per traitpair
    count <- 0
    for(i in traits) {
      for(j in traits) {
        traitpair <- paste(i,":",j,sep="",collapse=NULL)
        ij <- match(traitpair,alltraitpairs)
        count <- count + 1
        ci95lo <- object$variance.components[components,ij] - 1.96 * object$variance.components.se[components,ij]
        ci95hi <- object$variance.components[components,ij] + 1.96 * object$variance.components.se[components,ij]
        ctable <- data.frame(Traitpair=alltraitpairs[ij],
                     Estimate=object$variance.components[components,ij],
                     StdErr=object$variance.components.se[components,ij],
                     CI95lo=ci95lo,CI95hi=ci95hi,row.names=components)
        ctables[[count]] <- ctable
      }
    }
  }
  else {  # not bytrait

    ctables <- vector("list",length(components))  # one table per component
    count <- 0
    for(i in components) {
       count <- count + 1
       ci95lo <- object$variance.components[i,traitpairs] - 1.96 * object$variance.components.se[i,traitpairs]
       ci95hi <- object$variance.components[i,traitpairs] + 1.96 * object$variance.components.se[i,traitpairs]
       ctable <- data.frame(Component=i,
                    Estimate=object$variance.components[i,traitpairs],
                    StdErr=object$variance.components.se[i,traitpairs],
                    CI95lo=ci95lo,CI95hi=ci95hi,row.names=traitpairs)
       ctables[[count]] <- ctable
    }
  }
  retobj <- list(ctables=ctables,traits=traits, components=components, bytrait=bytrait, gls=gls, digits=digits)
  class(retobj) <- "csummarytables.dmm"
  return(retobj)

}
