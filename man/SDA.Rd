\name{SDA}
\alias{SDA}

\title{
Semi-parametric differential abuandance analysis
}
\description{
This function considers a two-part semi-parametric model for metabolomics and
proteomics data. A kernel-smoothed method is applied to estimate the regression
coefficients. And likelihood ratio test is constructed for differential
abundance analysis.
}
\usage{
SDA(sumExp)
}

\arguments{
  \item{sumExp}{
An object of 'SummarizedExperiment' class.
}

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
  \item{feat.names }{a vector of feature names }
}

\details{
    The differential abundance analysis is to compare metabolomic or proteomic
    profiles between different experimental groups, which utilizes a two-part
    model: a logistic regression model to characterize the zero proportion and a
    semi-parametric model to characterize non-zero values. Let
    \eqn{Y_{ig}}{Y_ig} be the random variable representing the abundance of
    feature \eqn{g} in subject \eqn{i}. This two-part model has the following
    form:

    \deqn{ \log(\frac{\pi_{ig}}{1-\pi_{ig}})=\gamma_{0g}+
    \boldsymbol{\gamma}_g \boldsymbol{X}_{i} }{ log(pi_ig/(1-pi_ig))=gamma_0g +
    gamma_g*X_i}
    \deqn{ \log(Y_{ig})=\boldsymbol{\beta}_g \boldsymbol{X}_i +
    \varepsilon_{ig} }{ log(Y_ig)=beta_g*X_i+ epsilon_ig }

    where \eqn{\pi_{ig}=Pr(Y_{ig}=0)}{pi_ig=Pr(Y_ig=0)} be the probability of
    point mass, \eqn{\boldsymbol{X}_i=(X_{i1},X_{i2},...,X_{iQ})^T}{X_i=(X_i1,
    X_i2,..., X_iQ)^T} is a Q-vector covariates that specifies the treatment
    conditions applied to subject \eqn{i}. The corresponding Q-vector of model
    parameters \eqn{\boldsymbol{\gamma}_g=(\gamma_{1g},\gamma_{2g},...,
    \gamma_{Qg})^T}{gamma_g=(gamma_1g, gamma_2g,...,gamma_Qg)^T} quantify the
    covariates effects on the fraction of zero values for feature \eqn{g} and
    \eqn{\gamma_{0g}}{gamma_0g} is the intercept. \eqn{\boldsymbol{\beta}_g=(
    \beta_{1g},\beta_{2g},...,\beta_{Qg})^T}{beta_g=(beta_1g, beta_2g,...,
    beta_Qg)
    ^T} is a Q-vector of model parameters quantifying the covariates effects on
    the non-zero values for the feature. And \eqn{\varepsilon_{ig}}{epsilon_ig}
    are independent error terms with a common but completely unspecified density
    function \eqn{f_g}{f_g}.

    Hypothesis testing on the effect of the \eqn{q}th covariate on the \eqn{g}th
    feature is performed by assessing \eqn{\gamma_{qg}}{gamma_qg} and \eqn{
    \beta_{qg}}{beta_qg}. Consider the null hypothesis \eqn{H_0}{H_0}: \eqn{
    \gamma_{qg}}{gamma_qg} and \eqn{\beta_{qg}}{beta_qg} against alternative
    hypothesis \eqn{H_1}{H_1}: at least one of the two parameters is non-zero.
    The p-value is calculated based on a chi-square distribution with 2 degrees
    of freedom. To adjust for multiple comparisons across features, the false
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
