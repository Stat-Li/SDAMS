\name{exampleSumExpset}
\alias{exampleSumExpset}
\docType{data}
\title{
A real 'SummarizedExperiment' object.
}
\description{
This example data is created based on a real proteomics data with 560 features
for 202 experimental subjects.
}
\usage{
  data(exampleSumExpset)
  }
\value{An object of \code{SummarizedExperiment} class.}

\seealso{
  \code{\link{SDA}}
}

\examples{
#------ load data --------
data(exampleSumExpset)
exampleSumExpset
feature = assay(exampleSumExpset) # access feature data
group = colData(exampleSumExpset)$grouping # access grouping information
SDA(exampleSumExpset)
}
\keyword{datasets}
