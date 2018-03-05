\name{createSEsetFromEnvir}
\alias{createSEsetFromEnvir}

\title{
SummarizedExperiment constructor
}
\description{
The constructor function \code{createMSsetFromEnvir} creates 
SummarizedExperiment object from R global environment.
}
\usage{
createSEsetFromEnvir(feature, group)
}

\arguments{
  \item{feature}{a matrix with row being features and column being subjects.}
  \item{group}{a phenotype data containing grouping information.}
}

\value{
An object of \code{SummarizedExperiment} class.
}

\author{
Yuntong Li <yuntong.li@uky.edu>,
Chi Wang <chi.wang@uky.edu>,
Li Chen <lichenuky@uky.edu>
}

\seealso{
  \code{\link{SDA}} input requires an object of \code{SummarizedExperiment} 
  class.
}

\examples{
set.seed(100)
featureInfo = matrix(runif(800,-2,5),ncol = 40)
featureInfo[featureInfo<0] = 0
rownames(featureInfo) = paste("feature",1:20,sep = '')
colnames(featureInfo) = paste('subject',1:40,sep = '')
groupInfo = data.frame(grouping=matrix(sample(0:1,40,replace = TRUE),ncol = 1))
rownames(groupInfo)=colnames(featureInfo)

exampleSEset = createSEsetFromEnvir(feature = featureInfo,group = groupInfo)
exampleSEset
}
