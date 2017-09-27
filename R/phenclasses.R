phenclasses <-
function(comcodes)
# phenclasses() - assemble all the phenotypic var classes for a setof comcodes
{
  if(length(comcodes) > 0) {
    nclasses <- length(comcodes[[1]])
    classes <- comcodes[[1]]
    if(length(comcodes) > 1) {
      for(kf in 2:length(comcodes) ) {
        nclasses <- nclasses * length(comcodes[[kf]])
        classes <- combpaste(classes,comcodes[[kf]])
      }
    }
  }
  return(classes)
}
