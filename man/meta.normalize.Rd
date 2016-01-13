\name{meta.normalize}
\alias{meta.Normalize}

\encoding{utf8}

\title{
immunoClust normalization step with the meta.clustering process 
}

\description{
Performs a normalization via linear regression of the
cell-cluster samples to the meta-clustering model. 
}

\usage{
meta.Normalize(P, N, K, W, M, S, G, Z, method=3)
}

\arguments{
\item{P}{The number of observed parameters for the cell event clusters.}
\item{N}{The number of cell-clustering experiments.}
\item{K}{The \eqn{N}-dimensional vector with the numbers of cell event
clusters in each experiment. The total number of clusters is
\eqn{totK = sum_{i=1}^K K_i}.}
\item{W}{The \eqn{totK}-dimensional vector with weights, i.e. number of 
events, of all clusters.}
\item{M}{The \eqn{totK x P}-dimensional matrix of all cluster means.}
\item{S}{The \eqn{totK x P x P}-dimensional matrix of all cluster covariance
matrices.}
\item{G}{The number of meta-clusters.}
\item{Z}{The \eqn{totK x G}-dimensional matrix with the A-Posterior 
probabilities for a cell-cluster belonging to a meta-cluster.}
\item{method}{
Alternative methods used for the normalization routine. Let \eqn{Y} denote 
the consensus meta-model build from all cell-event clusters of all experiments
using the A-posterior \eqn{Z} and \eqn{X} the cell-event clusters in each
experiment.

0   = no normalization

1   = \eqn{Y = a \times X}{ Y = a x X }

2   = \eqn{Y = a \times X + b}{ Y = a x X + b }

3   = \eqn{X = a \times Y}{ X = a x Y }

4   = \eqn{X = a \times Y + b}{ X = a x Y + b }

}
}

\details{
The regression used the cell-cluster and meta-cluster means weighted by the 
probabilities for a cell-cluster belonging to the meta-cluster. It builds
a consensus meta-model from all cell-clusters using the A-posterior 
probabilities \eqn{Z}. 
}

\value{
Returns the normalized cell-clusters means and co-variance matrices in a 
list-object with the following slots:
\describe{
\item{P}{The number of observed parameters for the cell event clusters.}
\item{N}{The number of cell-clustering experiments.}
\item{K}{The \eqn{N}-dimensional vector with the numbers of cell event
clusters in each experiment. The total number of clusters is
\eqn{totK = sum_{i=1}^K K_i}.}
\item{W}{The \eqn{totK}-dimensional vector with weights, i.e. number of 
events, of all clusters.}
\item{M}{The \eqn{totK x P}-dimensional matrix of all cluster means.}
\item{S}{The \eqn{totK x P x P}-dimensional matrix of all cluster covariance
matrices.}
}
}

\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}

\seealso{
\code{\link{meta.process}}, \code{\link{meta.Clustering}}
}

\examples{
data(dat.meta)
dat <- dat.meta$dat.clusters
res <- dat.meta$res.clusters
dat.norm <- meta.Normalize(dat$P, dat$N, dat$K, dat$clsEvents, 
    dat$M, dat$S, res@K, res@z)
}
