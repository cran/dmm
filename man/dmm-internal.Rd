\name{dmm-internal}
\alias{am.arel}
\alias{am.zandrel}
\alias{comtopar}
\alias{covcit}
\alias{covlyz}
\alias{covtopar}
\alias{covt}
\alias{crosseffectcovtopar}
\alias{dmesolve}
\alias{duplicated.first}
\alias{dyad.am.expect}
\alias{expect.v}
\alias{fixpaste}
\alias{gls.b.gmat}
\alias{gls.iter.b}
\alias{make.cohort}
\alias{partocov}
\alias{permpaste}
\alias{separ}
\alias{siga.posdef}
\alias{varlz}
\alias{varz}
\alias{vt}
\title{
Functions internal to dmm
}
\description{
Internal functions not documented and not visible to user.
}
\usage{
am.arel(df)
am.zandrel(mdf, df, k, l, v, x, y, cohortparts, components, relmat, ctable)
comtopar(v, l, siga, vara, vsiga, sesiga, ctable)
covcit(v, l, iv, vsiga, il, jl)
covlyz(cyz, y, z)
covtopar(covg, covt)
covt(v, l, vsiga, il, jl, kl)
crosseffectcovtopar(covx, covt, varcomp, var1, var2)
dmesolve(mdf, fixform = Ymat ~ 1, components = c("VarE(I)", "VarG(Ia)"),
     cohortform = NULL, posdef = T, gls = F,
     glsopt = list(maxiter = 200, bdamp = 0.8, stoptol = 0.01),
     dmeopt = "qr", ncomp.pcr = "rank", relmat = "inline", dmekeep = F, dmekeepfit = F)
duplicated.first(x)
dyad.am.expect(am, gls, dmeopt)
expect.v(am, siga, dme.explist)
fixpaste(x, fix)
gls.b.gmat(am, vinv, dme.explist, vmatblock)
gls.iter.b(am, start.b, start.siga, dyad.explist, glsopt, dmeopt, ctable,
           ncomp.pcr, dmekeepfit)
make.cohort(df, ce)
partocov(par, covt)
permpaste(x)
separ(siga, vsiga, v, l, fracta, correa)
siga.posdef(siga, am, ctable, varopt = "both", covopt = "bound")
varlz(vz, z)
varz(vlz, z)
vt(v, l, vsiga, il, jl)
}
%- maybe also 'usage' for other objects documented here.
\author{
Neville Jackson
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
