expect.v <-
function(am, siga, dme.explist){
# expect.v() - expected value of V matrix given siga estimates
    v <- matrix(0,am$l * am$n, am$l * am$n)
    for ( ic in 1 : nrow(siga)) {
      v <- v + kronecker(matrix(siga[ic, ],am$l,am$l),
                     matrix(dme.explist$vmat[ ,rownames(siga)[ic]],am$n,am$n),
                     make.dimnames=T)
    }

    return(v)
}
