\name{meta.regnorm}
\alias{meta.RegNorm}

\encoding{utf8}

\title{
immunoClust normalization procedure 
}

\description{
Performs a normalization via linear regression of the
sample clusters in x to the clusters in y. 
}

\usage{
meta.RegNorm(y, x, method=1, alpha=0.5)
}

\arguments{
\item{y}{ immunoClust-object with the destination clusters.}
\item{x}{ immunoClust-object woth the cluster to normalize.}
\item{method}{
Alternative methods used for the normalization routine. 

1   = \eqn{X = a \times Y}{ X = a x Y }

2   = \eqn{X = a \times Y + b}{ X = a x Y + b }
}
\item{alpha}{A value between 0 and 1 used to balance the bhattacharrya 
probabilities calculated with either the full covariance matrices or using only
the diagonal elements of it.}
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
\item{M}{The \eqn{totK x P}-dimensional matrix of all cluster means.}
\item{S}{The \eqn{totK x P x P}-dimensional matrix of all cluster covariance
matrices.}
}
}

\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}

\examples{
data(dat.meta)
data(dat.exp)
dat.norm <- meta.RegNorm(dat.meta$res.clusters, dat.exp[[1]])
}



