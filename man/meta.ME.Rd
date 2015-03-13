\name{meta.ME}
\alias{meta.ME}

\encoding{utf8}

\title{
immunoClust EM(t)-iteration on Cell-clusters
}

\description{
Performs an EM(t)-iteration on cell-clusters given an initial meta-cluster 
membership for the cell-clusters and returns the fitted meta-clusters 
information in an object of class 
\code{\link[=immunoClust.object]{immunoClust}}.
}

\usage{
meta.ME(P, N, K, W, M, S, label, B=500, tol=1e-5, method=20, bias=0.25, 
    alpha=1.0, min.class=0)
}

\arguments{
\item{P}{The number of observed parameters for the cell event clusters.}
\item{N}{The number of cell-clustering experiments.}
\item{K}{The \eqn{N}-dimensional vector with the numbers of cell event 
clusters in each experiment. The total number of clusters is 
\eqn{totK = sum_{i=1}^K K_i}.}
\item{W}{The \eqn{totK}-dimensional vector with weights, i.e. number of events,
of all clusters.}
\item{M}{The \eqn{totK x P}-dimensional matrix of all cluster means.}
\item{S}{The \eqn{totK x P x P}-dimensional matrix of all cluster covariance 
matrices.}
\item{label}{The \eqn{totK}-dimension integer vector with the initial 
cell-cluster to meta-cluster membership.}
\item{B}{The \eqn{totK x P x P}-dimensional matrix of all cluster covariance 
matrices.}
\item{tol}{The tolerance used to assess the convergence of the 
EM(t)-algorithms.}
\item{method}{
0   = KL-minimization not weighted

1   = BC-maximization not weighted

10  = BC-maximization weighted

2   = EMt-classification not weighted

20  = EMt-classification weighted
}
\item{bias}{The ICL-bias used in the EMt-iteration of the meta-clustering.}
\item{alpha}{A value between 0 and 1 used to balance the bhattacharrya 
probabilities calculated with either the full covariance matrices or using only
the diagonal elements of it. }
\item{min.class}{The minimum number of clusters for the final model.}
}

\details{
This function is used internally by the meta-clustering procedures 
\code{\link{meta.process}} and \code{\link{meta.Clustering}} in 
\emph{immuno}Clust.
}

\value{
The fitted meta-clusters information in an object of class 
\code{\link[=immunoClust.object]{immunoClust}}.
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
\code{\link{meta.process}}, \code{\link{meta.Clustering}}
}

\examples{
data(dat.exp)
d <- meta.exprs(dat.exp)
r <- meta.ME(d$P, d$N, d$K, d$clsEvents, d$M, d$S, label=rep(1,sum(d$K)))
}

\keyword{ cluster }
