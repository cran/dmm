\name{dmm-internal}
\alias{am.arel}
\alias{am.mcrel}
\alias{am.pyrel}
\alias{am.zandrel}
\alias{combpaste}
\alias{comtopar}
\alias{comtopar.specific}
\alias{countnotna}
\alias{covcit}
\alias{covlyz}
\alias{covtopar}
\alias{covt}
\alias{crossclasscovtopar}
\alias{crosseffectcovtopar}
\alias{csummary_specific}
\alias{dae.nonspecific.I}
\alias{dae.nonspecific}
\alias{dae.nonspecific.S}
\alias{dae.specific.I}
\alias{dae.specific}
\alias{dae.specific.S}
\alias{dmesolve}
\alias{dmm_array}
\alias{dmm_blockarray}
\alias{duplicated_first}
\alias{dyad.am.expect}
\alias{expect.v}
\alias{fixpaste}
\alias{genericvarvcnames}
\alias{genericvcnames}
\alias{gls.b.gmat}
\alias{gls.iter.b}
\alias{gsummary_specific}
\alias{is.cecov}
\alias{is.crossclass}
\alias{is.specific}
\alias{is.var}
\alias{is.withinclass}
\alias{make.cohort}
\alias{make.csummarytables}
\alias{make.gsummarytables}
\alias{make.matcline}
\alias{make.patyline}
\alias{match.cecov.nonspecific}
\alias{match.cecov.specific}
\alias{match.var.specific}
\alias{part.add}
\alias{partocov}
\alias{permpaste}
\alias{phenclasses}
\alias{selfpaste}
\alias{separ}
\alias{siga.posdef}
\alias{siga.posdef.specific}
\alias{sigatoie}
\alias{sigatovc}
\alias{traitpairstotraits}
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
am.mcrel(df)
am.pyrel(df)
am.zandrel(mdf, df, k, l, x, y, cohortparts, components,
     specific.components, relmat, ctable)
combpaste(x, y)
comtopar(v, l, siga, vara, vsiga, sesiga, ctable)
comtopar.specific(v, l, thisvc, vara, vthisvc, sethisvc, ctable,
     ic, phencovclass, longrownames, siga)
countnotna(x)
covcit(v, l, iv, vsiga, il, jl)
covlyz(cyz, y, z)
covtopar(covg, covt)
covt(v, l, vsiga, il, jl, kl)
crossclasscovtopar(covx, covt, var1, var2)
crosseffectcovtopar(covx, covt, varcomp, var1, var2)
csummary_specific(object, traitset, componentset, bytrait, gls, digits, ...)
dae.nonspecific.I(zpre, zpost, mmat, componentname, cnames, cnamesie,
     emat, vmat, icol, iecol, gls)
dae.nonspecific(zpre, rel, zpost, mmat, componentname, cnames, cnamesie,
     emat, vmat, icol, iecol, gls)
dae.nonspecific.S(zpre1, zpost1, zpre2, zpost2, zop, mmat,
     componentname, cnames, cnamesie, emat, vmat, icol, iecol, gls)
dae.specific.I(zpre, zpost, mmat, factorno, componentname,
     effnames, effcodes, effnandc, comcodes, varcodes, cnames, cnamesie,
     emat, vmat, icol, iecol, gls, ctable)
dae.specific(zpre, rel, zpost, mmat, factorno, componentname,
     effnames, effcodes, effnandc, comcodes, varcodes, cnames, cnamesie,
     emat, vmat, icol, iecol, gls, ctable)
dae.specific.S(zpre1, zpost1, zpre2, zpost2, zop, mmat, factorno,
     componentname, effnames, effcodes, effnandc, comcodes, varcodes,
     cnames, cnamesie, emat, vmat, icol, iecol, gls, ctable)
dmesolve(mdf, fixform = Ymat ~ 1, components = c("VarE(I)", "VarG(Ia)"),
     specific.components = NULL, cohortform = NULL, posdef = T, gls = F,
     glsopt = list(maxiter = 200, bdamp = 0.8, stoptol = 0.01),
     dmeopt = "qr", ncomp.pcr = "rank", relmat = "inline",
     dmekeep = F, dmekeepfit = F)
dmm_array(mdf, fixform, components, specific.components, cohortform,
     posdef, gls, glsopt, dmeopt, ncomp.pcr, relmat, dmekeep, dmekeepfit)
dmm_blockarray(mdf, fixform, components, specific.components, cohortform,
     posdef, gls, glsopt, dmeopt, ncomp.pcr, relmat, dmekeep, dmekeepfit, ...)
duplicated_first(x)
dyad.am.expect(am, gls, dmeopt)
expect.v(am, siga, dme.explist)
fixpaste(x, fix)
genericvarvcnames(x)
genericvcnames(x)
gls.b.gmat(am, vinv, dme.explist, vmatblock)
gls.iter.b(am, start.b, start.siga, dyad.explist, glsopt, dmeopt, ctable,
     ncomp.pcr, dmekeepfit)
gsummary_specific(dmmobj, traitset, componentset, bytrait, gls, digits, ...)
is.cecov(x, allcovs)
is.crossclass(x)
is.specific(x)
is.var(x)
is.withinclass(x)
make.cohort(df, ce)
make.csummarytables(object, traitset = "all",
     componentset = "all", bytrait = T, gls = F, digits = 3)
make.gsummarytables(dmmobj, traitset = "all",
     componentset = "all", bytrait = T, gls = F, digits = 3)
make.matcline(df)
make.patyline(df)
match.cecov.nonspecific(longrowname)
match.cecov.specific(longrowname)
match.var.specific(longrowname)
part.add(x, p, prow, pcol)
partocov(par, covt)
permpaste(x)
phenclasses(comcodes)
selfpaste(x)
separ(siga, vsiga, v, l, fracta, correa)
siga.posdef(siga, am, ctable, varopt = "both", covopt = "bound")
siga.posdef.specific(siga, am, ctable, varopt = "both", covopt = "bound")
sigatoie(cnames,cnamesie,siga,vsiga,sesiga,am,nsf)
sigatovc(siga, vsiga, sesiga, am, nsf)
traitpairstotraits(tp)
varlz(vz, z)
varz(vlz, z)
vt(v, l, vsiga, il, jl)
}
\author{
Neville Jackson
}
\keyword{ internal }
