\name{MSset-class}
\Rdversion{1.1}
\docType{class}
\alias{MSset-class}
\alias{MSset}

\alias{createMSsetFromEnvir}
\alias{createMSsetfromCSV}
\title{MSset object and constructors}
\description{
 The \code{MSset} is a s4 class used to store data from mass spectrometry. This
 class a subclass of \code{ExpressionSet}, with one more slot: group.
 The constructor functions \code{createMSsetFromEnvir} and
 \code{createMSsetfromCSV} create a \code{MSset} object from two kinds of
 situations: data in local environment or in csv files. See the vignette
 examples for each situation.
 }

\usage{
createMSsetFromEnvir(feature, group)

createMSsetfromCSV(featurePath, groupPath, rownames1=1,rownames2=1,
                  header1=TRUE,header2=TRUE)
}

\arguments{
  \item{feature}{a matrix with row being features and column being subjects.}

  \item{group}{a phenotype data containing grouping information.}

  \item{featurePath}{path for feature data.}

  \item{groupPath}{path for phenotype grouping data.}

  \item{rownames1}{indicator for feature data with row names. If NULL, row
  numbers are automatically generated.}

  \item{rownames2}{indicator for group data with row names. If NULL, row
  numbers are automatically generated.}

  \item{header1}{a logical value indicating whether the first row of feature is
  column names. The default value is TRUE.}

  \item{header2}{a logical value indicating whether the first row of group is
  column names. The default value is TRUE. If group input is a vector, set to
  False.}

}

\value{
     A \code{MSset} object.
}


\section{Methods}{
  \describe{
    \item{featuredata, featuredata<- :}{Access feature information.}
    \item{phenodata, phenodata<- :}{Access subject information.}
  }
}

\author{
    Yuntong Li <yuntong.li@uky.edu>,
    Chi Wang <chi.wang@uky.edu>,
    Li Chen <lichenuky@uky.edu>
}

\seealso{
  \code{\link{SDA}} input requires a class MSset.
}

\examples{
set.seed(100)
featureInfo = matrix(runif(800,-2,5),ncol = 40)
featureInfo[featureInfo<0]=0
rownames(featureInfo)=paste("feature",1:20,sep = '')
colnames(featureInfo)=paste('subject',1:40,sep = '')
groupInfo = data.frame(grouping=matrix(sample(0:1,40,replace = TRUE),ncol = 1))
rownames(groupInfo)=colnames(featureInfo)
exampleMSset2=createMSsetFromEnvir(feature = featureInfo,group = groupInfo)
exampleMSset2
head(featuredata(exampleMSset2)[,1:10])
head(phenotypedata(exampleMSset2))
}

\keyword{ class }

