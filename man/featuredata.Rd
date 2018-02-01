\name{featuredata}
\alias{featuredata,MSset-method}
\alias{featuredata<-,MSset,matrix-method}
\alias{featuredata}
\alias{featuredata<-}
\title{
  Acessor functions for the 'featuredata' slot in a MSset object.
}
\description{
  Basic feature information for metabolomics or proteomics data.

}
\usage{
\S4method{featuredata}{MSset}(object)
\S4method{featuredata}{MSset,matrix}(object) <- value
}

\arguments{

 \item{object}{A MSset object.}
 \item{value}{A matrix of feature data.}
}

\value{A matrix contains metabolomics or proteomics feature data.}


\author{
    Yuntong Li <yuntong.li@uky.edu>,
    Chi Wang <chi.wang@uky.edu>,
    Li Chen <lichenuky@uky.edu>
}

\seealso{
   \code{\link{phenotypedata}}
}
\examples{
data(exampleMSset)
# -------- access feature data ---------
fdata = featuredata(exampleMSset)
head(fdata)
# -------- assign a matrix -------------
featureInfo = matrix(runif(dim(fdata)[1]*dim(fdata)[2],-2,5),
                      ncol = dim(fdata)[2])
featureInfo[featureInfo<0]=0
featuredata(exampleMSset)=featureInfo
}

