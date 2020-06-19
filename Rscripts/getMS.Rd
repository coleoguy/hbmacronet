\name{getMS}
\alias{getMS}
\title{Creates an MS string}
\description{
This function takes a phylogeny and produces an MS ej string that can be used in other software.}
\usage{
getMS(tree, samps, report, n.site=NULL)
}
\arguments{
  \item{tree}{a phylo object}
  \item{samps}{number of coalescent simulations}
  \item{report}{type output for MS to produce. Currently limited to "T" for gene trees or "s" for segregating sites}
  \item{n.site}{This tells MS how many segregating sites to report - only used if "report=s".}
}
\details{
This function is designed to be used prior to the hybridization function.
}
\examples{
getMS(tree=rcoal(10),samps=2, report="T")
}


\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}
\keyword{ coalescent }
\keyword{ plotting phylogenies }
