\name{meta.SubClustering}
\alias{meta.SubClustering}
\alias{meta.TestSubCluster}

\encoding{utf8}

\title{
immunoClust Model Refinement Step in iterative Meta-clustering
}

\description{
These function tests each meta-cluster of a model for refining it into more 
sub-clusters and returns the refined cluster memberships in an integer array.s 
}

\usage{
meta.SubClustering(x, P, N, W, M, S, tol=1e-5, bias=0.25, thres=bias, 
    alpha=1.0, EM.method=20, verbose=FALSE)

meta.TestSubCluster(x, P, N, W, M, S, J=8, B=500, tol=1e-5, bias=0.5, 
    alpha=1.0, EM.method=2, HC.samples=2000) 
}

\arguments{
\item{x}{An immunoClust object with the initial model parameter 
    (\eqn{K}, \eqn{label}).}
\item{P}{The number of parameters.}
\item{N}{The number of clusters.}
\item{W}{The \eqn{N}-dimensional vector with cluster weights, i.e. numbers of 
    events in a cluster.}
\item{M}{The \eqn{N x P}-dimensional vector with cluster means.}
\item{S}{The \eqn{N x P x P}-dimensional vector with the cluster covariance 
    matrices.}
\item{tol}{The tolerance used to assess the convergence of the EM(t)-algorithms
    in Sub-Clustering.}
\item{bias}{he ICL-bias used in the EMt-algorithm.}
\item{thres}{Defines the threshold, below which an ICL-increase is meaningless.
    The threshold is given as the multiple (or fraction) of the costs for a 
    single cluster.}
\item{alpha}{A value between 0 and 1 used to balance the bhattacharrya 
    probabilities calculated with either the full covariance matrices or using 
    only the diagonal elements of it.}
\item{J}{The number of sub-models to be builded and tested for a particular 
    cluster.}
\item{B}{The maximum number of EM(t)-iterations in Sub-Clustering.}
\item{EM.method}{
0   = KL-minimization not weighted

1   = BC-maximization not weighted

10  = BC-maximization weighted

2   = EMt-classification not weighted

20  = EMt-classification weighted
}
\item{HC.samples}{The number of samples used for initial hierarchical 
    clustering.}
\item{verbose}{detailed messages during process}
}

\details{
These function are used internally by the meta-clustering procedures 
\code{\link{meta.process}} and \code{\link{meta.Clustering}} in 
\emph{immuno}Clust and are not intended to be used directly.
}

\value{
An integer array of length \eqn{N} containing the cell-clusters meta-cluster 
memberships of the refined model.
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
\code{\link{meta.process}}, \code{\link{meta.Clustering}}, 
\code{\link{meta.hclust}}
}

\examples{
data(dat.exp)
d <- meta.exprs(dat.exp)
#label <- rep(1,sum(d$K))
#label <- meta.SubClustering(d$P, sum(d$K), d$clsEvents, d$M, d$S, label=label)

r0 <- new("immunoClust", K=sum(d$K), label=rep(1,sum(d$K)))
label <- meta.SubClustering(r0, d$P, sum(d$K), d$clsEvents, d$M, d$S)

r1 <- meta.ME(d$P, d$N, d$K, d$clsEvents, d$M, d$S, label)


}

\keyword{ cluster }
