\name{SDAMS-package}
\docType{package}
\alias{SDAMS-package}

\title{
 SDAMS package for differential abundance analysis of Metabolomics, Proteomics 
 and single-cell RNA sequencing data
}
\description{
  SDAMS is an R package for differential abundance analysis of metabolomics and
  proteomics data from mass spectrometry and single-cell RNA sequencing data, 
  and the main function for differential abundance analysis
  is \code{\link{SDA}}. See the examples at \code{\link{SDA}} for
  basic analysis steps. SDAMS considers a two-part model, a logistic regression
  for the zero proportion and a semi-parametric log-linear model for the
  non-zero values.
}

\author{
  Yuntong Li <yuntong.li@uky.edu>,
  Chi Wang <chi.wang@uky.edu>,
  Li Chen <lichenuky@uky.edu>
}
\references{
Li, Y., Fan, T.W., Lane, A.N. et al. SDA: a semi-parametric differential
abundance analysis method for metabolomics and proteomics data.
BMC Bioinformatics 20, 501 (2019).
}
\keyword{package}
