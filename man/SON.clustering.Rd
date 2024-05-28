\name{meta.SON.clustering}
\alias{meta.SON.clustering}

\encoding{utf8}

\title{
meta clustering process with internal SON normalisation
}

\description{
The meta.SON.clustering is an extension of the meta-clustering
process co-clustering several samples cluster results. It integrates
a SON normalization step between the meta-clustering iterations.
}

\usage{
meta.SON.clustering(
meta,
cycles=6, alpha=0.5, scale.factor=2, scale.steps=0,
meta.iter=1, meta.bias=0.3, meta.thres=meta.bias, meta.tol=1e-5,
SON.cycles=1, SON.rlen=100, SON.deltas=c(1/SON.rlen,1/SON.rlen),
SON.blurring=c(2,0.1),
verbose=FALSE
)
}

\arguments{
\item{meta}{an immunoMeta-object for which the clustering should be
    refined.}
\item{cycles}{number of major iteration steps.}
\item{alpha}{The alpha value for calculation the bhattacharyya 
    probabilities.}
\item{scale.factor}{ scale factor for the internal model scaling step.}
\item{scale.steps}{ scale steps for the internal model scaling step. 0 means
    no model scaling.}
\item{meta.iter}{number of iterations for meta-clustering step }
\item{meta.bias}{ICL bias for meta-clustering step}
\item{meta.thres}{sub.thres for meta-clustering step}
\item{meta.tol}{maximal tolerance for meta-clustering step}
\item{SON.cycles}{number of cycles in SON normalization step.}
\item{SON.rlen}{runlength in SON normalization step}
\item{SON.deltas}{deltas parameter in SON normalization step}
\item{SON.blurring}{bluring parameter in SON normalisation step }
\item{verbose}{detailed messages during process}
}

\details{
For the refined \code{meta.SON.clustering} process a simple \code{meta.process}
should be performed first. The resulting \code{immunoMeta}-object then serves as
input data for the \code{meta.SON.clustering}.

Within the \code{meta.SON.clustering} between two meta.Clustering steps a 
SON normalization step is performed, which shifts the clusters of each sample 
towards the meta-clusters. The SON normalization for a sample consists of an 
optional first step to scale the model build by meta clusters best possible to 
the sample clusters. Afterwards, the meta clusters are moved to towards the 
sample clusters. This is done in a similar way to SOM clustering mapping. 
Finally, the sample clusters are retracted to the meta-clusters distribution. 
For this purpose the Bhattacharyya probabilities of sample and meta clusters
are used.
}
\value{
An immunoMeta-object for the co-clustering result.
}
\references{
pre-print
}
\author{
Till Sörensen \email{till.soerensen@bioretis.com}
}

\seealso{
\code{\link{meta.Clustering}}
}
\examples{
    data(dat.meta)
    meta <- meta.SON.clustering(dat.meta, cycles=2)
}
