\name{createSEsetFromCSV}
\alias{createSEsetFromCSV}

\title{
SummarizedExperiment constructor
}
\description{
The constructor function \code{createSEsetFromCSV} creates 
SummarizedExperiment object from csv.files.
}
\usage{
createSEsetFromCSV(featurePath, groupPath, rownames1=1,rownames2=1,
                  header1=TRUE,header2=TRUE)
}

\arguments{
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
directory1 <- system.file("extdata", package="SDAMS", mustWork=TRUE)
path1<-paste(directory1,"ProstateFeature.csv",sep="/")
directory2 <- system.file("extdata", package="SDAMS", mustWork=TRUE)
path2<-paste(directory2,"ProstateGroup.csv",sep="/")

exampleSEset = createSEsetFromCSV(path1,path2)
exampleSEset
}
