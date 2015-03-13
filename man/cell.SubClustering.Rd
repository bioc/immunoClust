\name{cell.SubClustering}
\alias{cell.SubClustering}
\alias{cell.TestSubCluster}

\encoding{utf8}

\title{
immunoClust Model Refinement Step in iterative Cell-events Clustering
}

\description{
These function tests each cell-cluster of a model for refining it into more 
sub-clusters and returns the refined model parameter in an object of class 
\code{\link[=immunoClust.object]{immunoClust}}.
}

\usage{
cell.SubClustering( x, dat, B=50, tol=1e-5, thres=0.1, bias=0.5,
                    sample.weights=1, sample.EM="MEt",
                    sample.number=1500, sample.standardize=TRUE,
                    extract.thres=0.8, modelName="mvt")

cell.TestSubCluster(x, y, t, cluster, J=8, B=500, tol=1e-5, bias=0.5,
                    sample.EM="MEt", sample.df=5, sample.number=1500, 
                    sample.standardize=TRUE, modelName="mvt") 
}

\arguments{
\item{x}{An immunoClust object with the initial model parameter 
(\eqn{K}, \eqn{w}, \eqn{mu}, \eqn{sigma}).}
\item{dat}{A numeric matrix, data frame of observations, or object of class 
flowFrame.} 
\item{B}{The maximum number of EM(t)-iterations in Sub-Clustering.}
\item{tol}{The tolerance used to assess the convergence of the EM(t)-algorithms
in Sub-Clustering.}
\item{thres}{Defines the threshold, below which an ICL-increase is meaningless.
The threshold is given as the multiple (or fraction) of the costs for a single 
cluster.}
\item{bias}{The ICL-bias used in the EMt-algorithm.}
\item{sample.weights}{Power of weights applied to hierarchical clustering, 
where the used weights are the probabilities of cluster membership.}
\item{sample.EM}{Used EM-algorithm; either \code{"MEt"} for EMt-iteration or 
\code{"ME"} for EM-iteration without test step.}
\item{sample.number}{The number of samples used for initial hierarchical 
clustering.}
\item{sample.standardize}{A numeric indicating whether the samples for 
hierarchical clustering are standardized (mean=0, SD=1).}
\item{extract.thres}{The threshold used for cluster data extraction.}
\item{modelName}{Used mixture model; either \code{mvt} for a \eqn{t}-mixture 
model or \code{mvn} for a Gaussian Mixture model.}

\item{y}{A numeric matrix of the observations beloning to the particular 
cluster. }
\item{t}{A numeric vector with the probability weights for the observations 
belonining to the particular cluster.}
\item{cluster}{An integer index of the particular cluster}
\item{J}{The number of sub-models to be builded and tested for a particular 
cluster. }
\item{sample.df}{Degree of freedom for the t-distibutions in a t-mixture model.
Has to be 5 in immunoClust.}
}

\details{
These function are used internally by the cell-clustering procedures of 
\code{\link{cell.process}} in \emph{immuno}Clust and are not intended to be used
directly. 
}

\value{
The cluster parameters of the refined model in an object of class 
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
\code{\link{cell.process}}, \code{\link{cell.hclust}}
}

\examples{
data(dat.fcs)
data(dat.exp)
dat.trans <- trans.ApplyToData(dat.exp[[1]], dat.fcs)
#need to re-calculate the cluster membership probabilities
# not stored in dat.exp
r1 <- cell.Classify(dat.exp[[1]], dat.trans)
summary(r1)
r2 <- cell.SubClustering(r1, dat.trans)
summary(r2)
}
\keyword{ cluster }
