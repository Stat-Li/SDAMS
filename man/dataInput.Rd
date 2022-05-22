\name{dataInput}
\alias{dataInput}

\alias{createSEFromCSV}
\alias{createSEFromMatrix}
\title{
    Mass spectrometry data input
}
\description{
    Two ways to input metabolomics or proteomics data from mass spectrometry or
    single-cell RNA sequencing data as \code{SummarizedExperiment}:

\enumerate{
\item \code{createSEFromCSV} creates SummarizedExperiment object from csv files;
\item \code{createSEFromMatrix} creates SummarizedExperiment object from
    separate matrices: one for feature/gene data and the other one for colData.
}

}
\usage{
createSEFromCSV(featurePath, colDataPath, rownames1 = 1, rownames2 = 1,
                  header1 = TRUE, header2 = TRUE)

createSEFromMatrix(feature, colData)
}

\arguments{
  \item{featurePath}{path for feature/gene data.}

  \item{colDataPath}{path for colData.}

  \item{rownames1}{indicator for feature/gene data with row names. If NULL, row
  numbers are automatically generated.}

  \item{rownames2}{indicator for colData with row names. If NULL, row
  numbers are automatically generated.}

  \item{header1}{a logical value indicating whether the first row of
  feature/gene is column names. The default value is TRUE.}

  \item{header2}{a logical value indicating whether the first row of colData is
  column names. The default value is TRUE. If colData input is a vector, set to
  False.}

  \item{feature}{a matrix with row being features/genes and column being
    subjects/cells.}

  \item{colData}{a column type data containing information about the
  subjects/cells.}
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
rownames(featureInfo) <- paste("gene", 1:20, sep = '')
colnames(featureInfo) <- paste('cell', 1:40, sep = '')
groupInfo <- data.frame(grouping=matrix(sample(0:1, 40, replace = TRUE),
                        ncol = 1))
rownames(groupInfo) <- colnames(featureInfo)

exampleSE <- createSEFromMatrix(feature = featureInfo, colData = groupInfo)
exampleSE
}
