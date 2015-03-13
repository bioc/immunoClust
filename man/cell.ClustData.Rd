\name{cell.ClustData}
\alias{cell.ClustData}

\encoding{utf8}

\title{
Model Based Clustering of Data for a pre-defined Number of Clusters
}
\description{
Performs EM-iteration on cell events, where an initial event cluster membership
is obtained by hierarchical clustering on a sample subset given a number of 
clusters.
}
\usage{
cell.ClustData(data, K, parameters=NULL, expName="immunoClust Experiment", 
                sample.seed=1, sample.number=1500, sample.standardize=TRUE,
                B=50, tol=1e-5, modelName="mvt") 
}
\arguments{
\item{data}{A numeric matrix, data frame of observations, or object of class 
flowFrame. Rows correspond to observations and columns correspond to measured 
parameters.}
\item{K}{Given number of clusters for the final model.}
\item{parameters}{A character vector specifying the parameters (columns) to be 
included in clustering. When it is left unspecified, all the parameters will be
used.}
\item{expName}{The name of the clustering experiment.}
\item{sample.seed}{The seed integer for the random number generator.}
\item{sample.number}{The maximum number of samples used for initial hierarchical
clustering.}
\item{sample.standardize}{A numeric indicating whether the samples for 
hierarchical clustering are standardized (mean=0, SD=1).}
\item{B}{The maximum number of EM-iterations.}
\item{tol}{The tolerance used to assess the convergence of the EM-algorithm.}
\item{modelName}{Used mixture model; either \code{"mvt"} for a t-mixture model 
or \code{"mvn"} for a Gaussian Mixture model.}
}
\details{
Although this function provides the possiblity to cluster an abitrary set of 
observed data into a fixed number of clusters, this function is used in the 
immunoClust-pipeline only for the calculation of the initial model with one 
cluster.
}
\value{
The fitted model cluster information in an object of class 
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
\code{\link[=immunoClust.object]{immunoClust-object}}, \code{\link{cell.hclust}}
}
\examples{
data(dat.fcs)
res <- cell.ClustData(dat.fcs, parameters=c("FSC-A", "SSC-A"), 5)
summary(res)
}
\keyword{ cluster }
