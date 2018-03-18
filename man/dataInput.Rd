\name{dataInput}
\alias{dataInput}

\alias{createSEFromCSV}
\alias{createSEFromMatrix}
\title{
    Mass spectrometry data input
}
\description{
    Two ways to input metabolomics or proteomics data from mass spectrometry as
\code{SummarizedExperiment}:

\enumerate{
\item \code{createSEFromCSV} creates SummarizedExperiment object from csv files;
\item \code{createSEFromMatrix} creates SummarizedExperiment object from
    separate matrices: one for feature data and the other one for group.
}

}
\usage{
createSEFromCSV(featurePath, groupPath, rownames1 = 1, rownames2 = 1,
                  header1 = TRUE, header2 = TRUE)

createSEFromMatrix(feature, group)
}

\arguments{
  \item{featurePath}{path for feature data.}

  \item{groupPath}{path for subject grouping data.}

  \item{rownames1}{indicator for feature data with row names. If NULL, row
  numbers are automatically generated.}

  \item{rownames2}{indicator for group data with row names. If NULL, row
  numbers are automatically generated.}

  \item{header1}{a logical value indicating whether the first row of feature is
  column names. The default value is TRUE.}

  \item{header2}{a logical value indicating whether the first row of group is
  column names. The default value is TRUE. If group input is a vector, set to
  False.}

  \item{feature}{a matrix with row being features and column being subjects.}

  \item{group}{a one column matrix containing subject grouping information.}
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
# ---------- csv input -------------
directory1 <- system.file("extdata", package = "SDAMS", mustWork = TRUE)
path1 <- file.path(directory1, "ProstateFeature.csv")
directory2 <- system.file("extdata", package = "SDAMS", mustWork = TRUE)
path2 <- file.path(directory2, "ProstateGroup.csv")

exampleSE <- createSEFromCSV(path1, path2)
exampleSE

# ---------- matrix input -------------
set.seed(100)
featureInfo <- matrix(runif(800, -2, 5), ncol = 40)
featureInfo[featureInfo<0] <- 0
rownames(featureInfo) <- paste("feature", 1:20, sep = '')
colnames(featureInfo) <- paste('subject', 1:40, sep = '')
groupInfo <- data.frame(grouping=matrix(sample(0:1, 40, replace = TRUE),
                        ncol = 1))
rownames(groupInfo) <- colnames(featureInfo)

exampleSE <- createSEFromMatrix(feature = featureInfo, group = groupInfo)
exampleSE
}
