\name{meta.clustering}
\alias{meta.Clustering}

\encoding{utf8}

\title{
Clustering of Cell-clusters in the immunoClust-pipeline
}

\description{
This function provides a direct access to the meta-clustering procedure. The 
method described and discussed in this manuscript is the EMt-classification 
(EM-method=20) with the number of events for each cluster as weights. It 
returns the fitted mixture model parameter in an object of class immunoClust.
}

\usage{
meta.Clustering(P, N, K, W, M, S, label=NULL, I.iter=10, B=500, tol=1e-5,
                bias=0.25, sub.thres = bias, alpha=0.5, EM.method=20,
                HC.samples=2000, norm.method=0, norm.blur=2, norm.minG=10,
                verbose=FALSE)
}

\arguments{
\item{P}{The number of observed parameters for the cell event clusters.}
\item{N}{The number of cell-clustering experiments.}
\item{K}{The \eqn{N}-dimensional vector with the numbers of cell event 
clusters in each experiment. The total number of clusters is  
\eqn{totK = sum_{i=1}^K K_i}.}
\item{W}{The \eqn{totK}-dimensional vector with weights of all clusters.}
\item{M}{The \eqn{totK x P}-dimensional matrix of all cluster means.}
\item{S}{The \eqn{totK x P x P}-dimensional matrix of all cluster covariance 
matrices.} 
\item{label}{Optional initial cluster assignment. If label equla NULL all
clusters are assigned in one cluster in the initial clustering step.}
\item{I.iter}{The maximum number of major iteration steps.}
\item{B}{The \eqn{totK x P x P}-dimensional matrix of all cluster covariance 
matrices.}
\item{tol}{The tolerance used to assess the convergence of the 
EM(t)-algorithms.}
\item{bias}{The ICL-bias used in the EMt-iteration of the meta-clustering.}
\item{sub.thres}{Defines the threshold, below which an ICL-increase is
meaningless. The threshold is given as the multiple (or fraction) of the costs
for a single cluster.}
\item{alpha}{A value between 0 and 1 used to balance the bhattacharrya 
probabilities calculated with either the full covariance matrices or using only
the diagonal elements of it. When working with uncompensated FC data very high 
correlations between parameters may be observed due to spill over. This leads 
to a very low bhattacharrya probability for two clusters even if they are 
located nearby. Using a mixture of the probabilities calculated with the 
complete covariance matrices and the variance information of each parameter 
avoids this problem. With a value of alpha=1, only the probabilities with 
complete covariance matrices are applied. A reasonable value for alpha is 0.5.}
\item{EM.method}{
0   = KL-minimization not weighted

1   = BC-maximization not weighted

10  = BC-maximization weighted

2   = EMt-classification not weighted

20  = EMt-classification weighted
}
\item{HC.samples}{The number of samples used for initial hierarchical
    clustering.}
\item{norm.method}{Normalization function; see \code{\link{meta.Normalize}} for
details.}
\item{norm.blur}{For the normalization step the a-posteriori probabilites of the
cell-clusters belonging to a meta.clusters a used. In order to capture narrow 
cell-clusters reasonable the co-variance of the cell-clusters is blured for the
a-posteriori probabilities in the normalization step.}
\item{norm.minG}{Minimum number of obtained meta-clusters required to process 
the normalization step in the major iteration loop.}
\item{verbose}{detailed messages during process}
}

\details{
This function is used internally by the meta-clustering procedure 
\code{\link{meta.process}} in \emph{immuno}Clust.
}

\value{
The fitted model information in an object of class 
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
\code{\link[=immunoClust.object]{immunoClust-object}}, 
\code{\link{meta.SubClustering}}, \code{\link{meta.process}}
}

\examples{
data(dat.exp)
d <- meta.exprs(dat.exp)
res <- meta.Clustering(d$P, d$N, d$K, d$clsEvents, d$M, d$S)
}

\keyword{ cluster }
