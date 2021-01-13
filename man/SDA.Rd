\name{SDA}
\alias{SDA}

\title{
Semi-parametric differential abuandance/expression analysis
}
\description{
This function considers a two-part semi-parametric model for metabolomics,
proteomics and single-cell RNA sequencing data. A kernel-smoothed method is
applied to estimate the regression
coefficients. And likelihood ratio test is constructed for differential
abundance/expression analysis.
}
\usage{
SDA(sumExp, VOI = NULL, ...)
}

\arguments{
  \item{sumExp}{An object of 'SummarizedExperiment' class.}
  \item{VOI}{Variable of interest. Default is NULL, when there is only one
  covariate, otherwise it must be one of the column names in colData.}
  \item{...}{Additional arguments passed to \code{\link{qvalue}}.}

}

\value{
A list containing the following components:
  \item{gamma }{a vector of point estimators for \eqn{\gamma_g}{gamma_g} in the
  logistic model (binary part) }
  \item{beta }{a vector of point estimators for \eqn{\beta_g}{beta_g} in the
  semi-parametric model (non-zero part) }
  \item{pv_gamma }{a vector of one-part p-values for \eqn{\gamma_g}{gamma_g} }
  \item{pv_beta }{a vector of one-part p-values for \eqn{\beta_g}{beta_g} }
  \item{qv_gamma }{a vector of one-part q-values for \eqn{\gamma_g}{gamma_g} }
  \item{qv_beta }{a vector of one-part q-values for \eqn{\beta_g}{beta_g} }
  \item{pv_2part }{a vector of two-part p-values for overall test }
  \item{qv_2part }{a vector of two-part q-values for overall test }
  \item{feat.names }{a vector of feature/gene names }
}

\details{
    The differential abundance/expression analysis is to compare metabolomic or             proteomic profiles or gene expression between different experimental groups,
    which utilizes a two-part
    model: a logistic regression model to characterize the zero proportion and a
    semi-parametric model to characterize non-zero values. Let
    \eqn{Y_{i}}{Y_i} be the random variable and \eqn{X_{i}}{X_i} is a vector of             covariates. This two-part model has the following
    form:

    \deqn{ \log(\frac{\pi_{i}}{1-\pi_{i}})=\gamma_{0}+
    \boldsymbol{\gamma} \boldsymbol{X}_{i} }{ log(pi_i/(1-pi_i))=gamma_0 +
    gamma*X_i}
    \deqn{ \log(Y_{i})=\boldsymbol{\beta} \boldsymbol{X}_i +
    \varepsilon_{i} }{ log(Y_i)=beta*X_i+ epsilon_i }

    where \eqn{\pi_{i}=Pr(Y_{i}=0)}{pi_i=Pr(Y_i=0)}. The model parameters
    \eqn{\boldsymbol{\gamma}}{gamma} quantify the covariates effects on the
    fraction of zero values and \eqn{\gamma_0}{gamma_0} is the intercept.
    \eqn{\boldsymbol{\beta}}{beta} are the model parameters quantifying the
    covariates effects on the non-zero values, \eqn{\varepsilon_{i}}{epsilon_i}
    are independent error terms with a common but completely unspecified density            function \eqn{f}{f}.

    For differential abundant analysis on data from mass spectrometry,
    \eqn{Y_{i}}{Y_i} represents the abundance of certain feature for subject
    \eqn{i}{i}, \eqn{\pi_{i}}{pi_i} is the probability of point mass.
    \eqn{\boldsymbol{X}_i=(X_{i1},X_{i2},...,X_{iQ})^T}{X_i=(X_i1,
    X_i2,..., X_iQ)^T} is a Q-vector of covariates that specifies the treatment
    conditions applied to subject \eqn{i}. The corresponding Q-vector of model
    parameters \eqn{\boldsymbol{\gamma}=(\gamma_{1},\gamma_{2},...,
    \gamma_{Q})^T}{gamma=(gamma_1, gamma_2,...,gamma_Q)^T} and
    \eqn{\boldsymbol{\beta}=(
    \beta_{1},\beta_{2},...,\beta_{Q})^T}{beta=(beta_1, beta_2,...,
    beta_Q)^T} quantify the covariates effects for certain feature. Hypothesis
    testing on the effect of the \eqn{q}th covariate on certain
    feature is performed by assessing \eqn{\gamma_{q}}{gamma_q} and \eqn{
    \beta_{q}}{beta_q}. Consider the null hypothesis \eqn{H_0}{H_0}: \eqn{
    \gamma_{q}=0}{gamma_q=0} and \eqn{\beta_{q}=0}{beta_q=0} against alternative
    hypothesis \eqn{H_1}{H_1}: at least one of the two parameters is non-zero.
    We also consider the hypotheses for testing \eqn{\gamma_{q}=0}{gamma_q=0} and
    \eqn{\beta_{q}=0}{beta_q=0} individually.


    For differential expression analysis on single-cell RNA sequencing data,
    \eqn{Y_{i}}{Y_i} represents represents the expression (TPM value) of certain gene in     \eqn{i}{i}th cell, \eqn{\pi_{i}}{pi_i} is the drop-out probability.
    \eqn{\boldsymbol{X}_i=(Z_{i}, \boldsymbol{W}_i)^T}{X_i=(Z_i,
    W_i)^T} is a vector of covariates with \eqn{Z_i}{Z_i} being a binary
    indicator of the cell population under comparison and
    \eqn{\boldsymbol{W}_i}{W_i} being a vector of other covariates, e.g. cell
    size, and \eqn{\boldsymbol{\gamma}=(\gamma_{Z},\boldsymbol{\gamma}_W)}{gamma            =(gamma_Z, gamma_W)} and
    \eqn{\boldsymbol{\beta}=(\beta_{Z},\boldsymbol{\beta}_W)}{beta=
    (beta_Z, beta_W)} are model parameters. Hypothesis
    testing on the effect of different cell subpopulations on certain
    gene is performed by assessing \eqn{\gamma_{Z}}{gamma_Z} and \eqn{
    \beta_{Z}}{beta_Z}. For
    each gene, the likelihood ratio test is performed on the null hypothesis
    \eqn{H_0}{H_0}: \eqn{
    \gamma_{Z}=0}{gamma_Z=0} and \eqn{\beta_{Z}=0}{beta_Z=0} against alternative
    hypothesis \eqn{H_1}{H_1}: at least one of the two parameters is non-zero.
    We also consider the hypotheses for testing \eqn{\gamma_{Z}=0}{gamma_Z=0} and
    \eqn{\beta_{Z}=0}{beta_Z=0} individually.

    The p-value is calculated based on an asympotic chi-squared distribution. To
    adjust for multiple comparisons across features, the false
    discovery discovery rate (FDR) q-value is calculated based on the
    \link[qvalue]{qvalue} function in R/Bioconductor.

}

\author{
Yuntong Li <yuntong.li@uky.edu>,
Chi Wang <chi.wang@uky.edu>,
Li Chen <lichenuky@uky.edu>
}



\examples{
##--------- load data ------------
data(exampleSumExp)

results = SDA(exampleSumExp)

##------ two part q-values -------
results$qv_2part
}

\keyword{ model }
