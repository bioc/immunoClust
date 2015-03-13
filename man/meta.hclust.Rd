\name{meta.hclust}
\alias{meta.hclust}

\encoding{utf8}

\title{
Hierarchical Meta-clustering of Cell-clusters in the immunoClust-pipeline
}

\description{
Performs agglomerative clustering on cell-clusters. It is used in the interative
meta-clustering approach of \emph{immuno}Clust to obtain an initial meta-cluster
membership for the EM(t)-iteration.
}

\usage{
meta.hclust(P, N, W, M, S)
}

\arguments{
\item{P}{The number of parameters.}
\item{N}{The number of clusters.}
\item{W}{The \eqn{N}-dimensional vector with cluster weights, i.e. numbers of 
events in a cluster.}
\item{M}{The \eqn{N \times P}{N x P}-dimensional vector with cluster means.}
\item{S}{The \eqn{N \times P \times P}{N x P x P}-dimensional vector with 
cluster covariance matrices.}
}

\details{
This function is used internally in \code{\link{meta.TestSubCluster}} of 
\pkg{immunoClust}.
}

\value{
A numeric \eqn{(N-1) \times 2}{(N-1) x 2}-dimensional matrix which gives the 
minimum index for observations in each of the two clusters merged at the 
\emph{i}th step in each row.
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

\note{
The merging distances need not to be monotonic increasing.
}

\seealso{
\code{\link{meta.TestSubCluster}}, \code{\link{meta.process}}
}

\examples{
data(dat.exp)
r <- dat.exp[[1]]
hcPairs <- meta.hclust(r@P, r@K, r@w, r@mu, t(apply(r@sigma,1,c)))
}

\keyword{ cluster }
