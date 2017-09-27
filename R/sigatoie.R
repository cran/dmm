sigatoie <-
function(cnames,cnamesie,siga,vsiga,sesiga,am,nsf)
#  sigatoie() - add inestimable rows to siga,sesiga
#             - add inestimable rows and cols to vsiga
{
  fullv <- length(cnames) + length(cnamesie)
  sigaie <- matrix(NA,fullv,ncol(siga))
  sesigaie <- matrix(NA,fullv,ncol(siga))
  vsigaie <- matrix(NA,fullv*ncol(siga),fullv*ncol(siga))
  rownames(sigaie) <- c(cnames,cnamesie)
  colnames(sigaie) <- colnames(siga)
  rownames(sesigaie) <- c(cnames,cnamesie)
  colnames(sesigaie) <- colnames(sesiga)
  rownames(vsigaie) <-  combpaste(colnames(siga),c(cnames,cnamesie))
  colnames(vsigaie) <- rownames(vsigaie)
  sigaie[cnames,] <- siga
  sesigaie[cnames,] <- sesiga
  vsigaie[rownames(vsiga),colnames(vsiga)] <- vsiga
  outlist <- list(sigaie =sigaie ,sesigaie=sesigaie,vsigaie=vsigaie)
  return(outlist)
}
