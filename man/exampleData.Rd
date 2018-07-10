\name{exampleData}
\alias{exampleData}
\alias{exampleSumExp}
\docType{data}
\title{
    An example data for the SDAMS package
}
\description{
SDAMS package provides two different data formats for prostate cancer proteomics
data, which is from the human urinary proteome database(\url{
http://mosaiques-diagnostics.de/mosaiques-diagnostics/human-urinary-proteom-
database}). There are 526 prostate cancer subjects and 1503 healthy subjects. A
total of 5605 proteomic features were measured for each subject. For
illustration purpose, we took a 10\% subsample randomly from this real data.
This example data contains 560 proteomic features for 202 experimental subjects
with 49 prostate cancer subjects and 153 healthy subjects.
\code{exampleSumExp.rda} is an object of \code{SummarizedExperiment} class which
stores the information of both proteomic features and experimental subjects.
\code{ProstateFeature.csv} contains a matrix-like proteomic feature data and
\code{ProstateGroup.csv} contains a single column of experimental subject group
data.
}
\usage{
  data(exampleSumExp)
  }
\value{An object of \code{SummarizedExperiment} class.}

\seealso{
  \code{\link{SDA}}
}

\examples{
#------ load data --------
data(exampleSumExp)
exampleSumExp
feature = assay(exampleSumExp) # access feature data
group = colData(exampleSumExp)$grouping # access grouping information
SDA(exampleSumExp)

}
\references{
Siwy, J., Mullen, W., Golovko, I., Franke, J., and Zurbig, P. (2011). Human
urinary peptide database for multiple disease biomarker discovery.
PROTEOMICS-Clinical Applications 5, 367-374.
}
\keyword{datasets}
