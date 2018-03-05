\name{SDA}
\alias{SDA}

\title{
semi-parametric differential abuandance analysis
}
\description{
This function considers a two-part semi-parametric model for metabolomics and
proteomics data. A kernel-smoothed method is applied to estimate the regression
coefficients. And likelihood ratio test is constructed for differential
abundance analysis.
}
\usage{
SDA(sumExpset)
}

\arguments{
  \item{sumExpset}{
An object of 'SummarizedExperiment' class.
}

}

\value{
A list containing the following components:
  \item{gamma }{a vector of point estimators for binary part }
  \item{beta }{a vector of ponit estimators for non-zero part }
  \item{pv_gamma }{a vector of p-values for gamma }
  \item{pv_beta }{a vector of p-values for beta }
  \item{qv_gamma }{a vector of q-values for gamma }
  \item{qv_beta }{a vector of q-values for beta }
  \item{pv_2part }{a vector of p-values for the combined two-part }
  \item{qv_2part }{a vector of q-values for the combined two-part }
  \item{feat.names }{a vector of feature names }
}

\author{
Yuntong Li <yuntong.li@uky.edu>,
Chi Wang <chi.wang@uky.edu>,
Li Chen <lichenuky@uky.edu>
}



\examples{
##--------- load data ------------
data(exampleSumExpset)

results = SDA(exampleSumExpset)
results$qv_2part
}

\keyword{ model }
