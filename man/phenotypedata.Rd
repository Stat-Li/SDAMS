\name{phenotypedata}
\alias{phenotypedata,MSset-method}
\alias{phenotypedata}
\alias{phenotypedata<-}
\alias{phenotypedata<-,MSset,data.frame-method}
\title{
  Acessor functions for the 'phenoData' slot in a MSset object.
}
\description{
  Basic phenotype information for metabolomics or proteomics data.

}
\usage{
\S4method{phenotypedata}{MSset}(object)
\S4method{phenotypedata}{MSset,data.frame}(object) <- value
}

\arguments{

 \item{object}{A MSset object.}
 \item{value}{A phenotype data.}
}

\value{A one column matrix contains metabolomics or proteomics subject grouping
      data.}


\author{
    Yuntong Li <yuntong.li@uky.edu>,
    Chi Wang <chi.wang@uky.edu>,
    Li Chen <lichenuky@uky.edu>
}

\seealso{
   \code{\link{featuredata}}
}
\examples{
data(exampleMSset)
# ----------- access subject grouping data ------------
pdata = phenotypedata(exampleMSset)
pdata
# -------- assign a one column data frame -------------
groupInfo = data.frame(grouping=matrix(sample(0:1,dim(pdata)[1],replace = TRUE),
                        ncol = 1))
phenotypedata(exampleMSset)=groupInfo
}

