\name{meta.Normalize}
\alias{meta.Normalize}
\alias{meta.GPA}
\alias{meta.Scale}

\encoding{utf8}

\title{
immunoClust normalization step with the meta.clustering process 
}

\description{
Performs a normalization vio linear or quadratic regression of the
cell-cluster samples to the meta-clustering model. The meta.GPA function
is in an experimental state and processes procrustes analysis instead of
linear regression.
}

\usage{
meta.Scale(P, N, K, W, M, S, method=5)

meta.GPA(P, N, K, W, M, S, G, Z)

meta.Normalize(P, N, K, W, M, S, G, Z, degree=1)
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
\item{Z}{The \eqn{totK x G}-dimensiona matrix with the A-Posterior 
probabilities for a cell-cluster belonging to a meta-cluster.}
\item{method}{Experimental alternative methods used for the scaling routine:

0   = sample median absolute deviation to average median absolute deviation

1   = sample trimed sd (mean=0) to average trimed sd (mean=0)

2   = sample trimed mean/sd to average trimed mean/sd

3   = sample trimed mean/sd to mean=0, sd=1

4   = sample weighted mean/sd to average weighted mean/sd

5   = 0.9-quantile of sample means to average 0.9-quantiles of sample means
}
\item{degree}{The degree of the regression ploynom used for the normalization.}
}

\details{
The normalization performs a weighted regression through the origin. Thus, 
only a scaling without a translation for the cell-clusters in one sample is 
applied.
The regression used the cell-cluster and meta-cluster means weighted by the 
probabilities for a cell-cluster belonging to the meta-cluster.

The meta.Scale function is used only for the first sample cell-clusters 
scaling within the normalization iteration in the meta.clustering process.

The meta.GPA function is experimental and not used for normalization jet.
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
