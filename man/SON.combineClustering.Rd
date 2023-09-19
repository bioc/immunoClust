\name{meta.SON.combineClustering}
\alias{meta.SON.combineClustering}

\encoding{utf8}

\title{
Transfer the annotation of an immunoMeta-object to an
immunoClust-object.
}
\description{
An immunoMeta-object is co-clustered with an immunoClust-object 
of the same parameter structure. Co-clustering includes SON normalization 
steps. The returned immnuoCLust-object contians the meta-clusters unchanged 
in order and numeration. 
}

\usage{
meta.SON.combineClustering(
meta, res, par=seq_len(npar(meta)),
map.cluster=seq_len(ncls(meta)),
use.cluster=seq_len(ncls(res)),
meta.alpha=0.5, meta.bias=0.1, meta.iter=100, meta.tol=1e-5,
SON.method=1, SON.cycles=4, SON.rlen=10,
SON.deltas=c(1/SON.rlen,1/SON.rlen), SON.blurring=c(2,1),
traceG=c(), traceK=c())
}

\arguments{
\item{meta}{The annotated immunoMeta-object.}
\item{res}{An immunoClust-object as results from cell-event
    clustering for a sample}
\item{par}{An integer array with the parameters to be used for
    SON mapping.}
\item{map.cluster}{The model clusters to be used for SON mapping.}
\item{use.cluster}{the sample clusters to be used for SON mapping.}
\item{meta.alpha}{The alpha value in calculation the bhattacharyya 
    probabilities.}

\item{meta.bias}{The ICL bias for meta co-clustering step.}
\item{meta.iter}{Maximal iterations in the meta co-clustering step.}
\item{meta.tol}{Maximal tolerance for meta co-clustering step.}
\item{SON.method}{Method selection for SON normalization step.}
\item{SON.cycles}{Number cycles in the SON normalization step. }
\item{SON.rlen}{ runlength in the SON normalization step.}
\item{SON.deltas}{ delta parameter in the SON normalization step.}
\item{SON.blurring}{ blurring parameter in the SON normalization step.}
\item{traceG}{An array of model cluster to trace in the process.}
\item{traceK}{An array of sample cluster to trace in the process. }

}
\details{
The co-clustering consists of a normalization and meta-clustering step.
A sample cluster is than labeled according to its corresponding meta
cluster. The SON-normalization and meta-clustering steps are parameterised by 
the SON and meta arguments.
}
\value{
An immunoClust-object from meta-clusters and combined
observation from meta- and samples-cluster. The first
G elements of the label coresponds to the meta-clusters, afterwards
the labelling of the samples-clusters indicates the nearest 
meta-cluster for the sample-cluster.
}
\references{
in progress
}
\author{
Till Sörensen \email{till.soerensen@bioretis.com}
}

\seealso{
\code{\link{meta.Clustering}}
}
\examples{
    data(dat.exp)
    data(dat.meta)
    res <- meta.SON.combineClustering(dat.meta, dat.exp[[1]], SON.cycles=2)
}
