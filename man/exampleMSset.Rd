\name{exampleMSset}
\alias{exampleMSset}
\docType{data}
\title{
A real 'MSset' object.
}
\description{
This example data is created based on a real proteomics data with 560 features
for 202 experimental subjects.
}
\usage{
  data(exampleMSset)
  }
\value{An object of MSset.}

\seealso{
  \code{\link{SDA}}
}

\examples{
#------ load data --------
data(exampleMSset)
exampleMSset
SDA(exampleMSset)
}
\keyword{datasets}
