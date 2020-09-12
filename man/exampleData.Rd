\name{exampleData}
\alias{exampleData}
\alias{exampleSumExp}
\alias{exampleSingleCell}
\docType{data}
\title{
    Two example datasets for SDAMS package
}
\description{
SDAMS package provides two types of example datasets: one is prostate cancer
proteomics data and the other one is single cell RNA sequencing data.
\enumerate{
\item  For prostate cancer proteomics data, it is from the human urinary
proteome database(\url{
http://mosaiques-diagnostics.de/mosaiques-diagnostics/human-urinary-proteom-
database}). There are 526 prostate cancer subjects and 1503 healthy subjects. A
total of 5605 proteomic features were measured for each subject. For
illustration purpose, we took a 10\% subsample randomly from this real data.
This example data contains 560 proteomic features for 202 experimental subjects
with 49 prostate cancer subjects and 153 healthy subjects. SDAMS package
provides two different kinds of data formats for prostate cancer proteomics
data.
\code{exampleSumExp.rda} is an object of \code{SummarizedExperiment} class which
stores the information of both proteomic features and experimental subjects.
\code{ProstateFeature.csv} contains a matrix-like proteomic feature data and
\code{ProstateGroup.csv} contains a single column of experimental subject group
data.
\item For single cell RNA sequencing data, it is in the form of transcripts per
kilobase million (TPM). The count data can be found at Gene Expression Omnibus
(GEO) database with Accession No. GSE29087. There are 92 single
cells (48 mouse embryonic stem (ES) cells and 44 mouse embryonic fibroblasts
(MEF)) that were analyzed. The example data provided by SDAMS contains 10\% of
genes which are randomly sampled from the raw dataset.
\code{exampleSingleCell.rda} is an object of \code{SummarizedExperiment} class
which stores the information of both gene expression and cell information.
}

}
\usage{
  data(exampleSumExp)
  data(exampleSingleCell)
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

Islam, S., Kjallquist, U., Moliner, A., Zajac, P., Fan, J. B., Lonnerberg, P., & Linnarsson, S. (2011). Characterization of the single-cell transcriptional
landscape by highly multiplex RNA-seq. Genome research, 21(7), 1160-1167.
}
\keyword{datasets}
