\name{meta.process}
\alias{meta.process}

\encoding{utf8}

\title{
Meta-clustering of Cell-clusters in the immunoClust-pipeline
}

\description{
This function performs iterative model based clustering on the clusters obtained
by \code{\link{cell.process}} of several samples. Its input is a vector of the 
\code{immunoClust-objects} of the samples.
}

\usage{
meta.process(exp, dat.subset=c(), meta.iter=10, tol=1e-05, meta.bias=0.2,
            meta.alpha=.5, norm.method=0, norm.blur=2, norm.minG=10)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{exp}{A vector of \code{list} objects, each \code{list} contains the 
cell-clustering result of a sample in the \code{res} field. Addition fields are
\code{name} and \code{fsc} containing the cell-sample name and fcs-filename, 
which are used for data output and plot routines.}
\item{dat.subset}{A numeric vector defining the used observed parameters for the
meta-clustering. If unset, all parameters in the cell-clustering results are 
used.}
\item{meta.iter}{The number of major iterations.}
\item{tol}{The tolerance used to assess the convergence of the 
EM(t)-algorithms.}
\item{meta.bias}{The ICL-bias used in the EMt-iteration of the meta-clustering.}
\item{meta.alpha}{A value between 0 and 1 used to balance the bhattacharrya 
probabilities calculated with either the full covariance matrices or using only
the diagonal elements of it. When working with uncompensated FC data, very high
correlations between parameters may be observed due to spill over. This leads to
a very low bhattacharrya probability for two clusters even if they are located 
nearby. Using a mixture of the probabilities calculated with the complete 
covariance matrices and the variance information of each parameter avoids this 
problem. With a value of alpha=1, only the probabilities with complete 
covariance matrices are applied. A reasonable value for alpha is 0.5.}
\item{norm.method}{A numeric selector for the normalization step to be
performed during the major iteration.} 
\item{norm.blur}{The bluring constant by which the cell-clusters co-variance 
matrices are increased within the normalization step.} 
\item{norm.minG}{Minimum number of meta-clusters required before processing
the normalization step.}

}

\value{
The function returns a \code{\link[=immunoMeta.object]{immunoMeta}} with the 
following components:
\tabular{ll}{
\code{dat.clusters} \tab A \code{dat} list-object of the cell event clusters 
used for meta-clustering.\cr
\code{res.clusters} \tab The 
\code{\link[=immunoClust.object]{immunoClust-object}} of the fitted 
meta-clustering mixture model.\cr
\code{dat.scatter}  \tab A \code{dat} list-object of the scatter parameters for
the cell event clusters used for scatter clustering.\cr
\code{res.scatter}  \tab The 
\code{\link[=immunoClust.object]{immunoClust-object}} of the fitted 
scatter-clustering mixture model.\cr
\code{gating}       \tab A list-object containing the hierarchical 
gating-tree.\cr
}


}

\references{
Sörensen, T., Baumgart, S., Durek, P., Grützkau, A. and Häupl, T. 
immunoClust - an automated analysis pipeline for the identification of 
immunophenotypic signatures in high-dimensional cytometric datasets. 
\emph{Cytometry A} (accepted). 
}
\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}

\seealso{
\code{\link[=immunoMeta.object]{immunoMeta-object}},
\code{\link[=immunoClust.object]{immunoClust-object}}, 
\code{\link{meta.Clustering}}, \code{\link{meta.export}}, 
\code{\link{cell.process}}
}

\examples{
data(dat.exp)
meta <- meta.process(dat.exp)
summary(meta)
tbl <- meta.numEvents(meta)
}

\keyword{ cluster }
