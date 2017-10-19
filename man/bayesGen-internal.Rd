% --- bayesGen-internal.Rd ---

\name{bayesGen-internal}
\alias{bayesGen-internal}
\alias{getInfoBayesSNPassoc}
\alias{getInfoBayesCNVassoc}
\alias{snpMatrix2bayes}
\alias{datChr}
\alias{writeFileGenoDat}

\title{Internal bayesGen functions}
\description{Internal bayesGen functions}
\usage{
getInfoBayesSNPassoc(x, Ngroups, Nvar, names.groups, names.SNPs)
snpMatrix2bayes(x)
datChr(x, geno, annotation)
getInfoBayesCNVassoc(x, N.groups, N.cnvs, names.groups, names.CNVs)
writeFileGenoDat(gg, allesum=FALSE, caco)
}

\details{These are not to be called by the user}
\keyword{internal}
