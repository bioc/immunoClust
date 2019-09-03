\name{immunoMeta-class}
\docType{class}
\alias{immunoMeta-class}
\alias{immunoMeta-object}
\alias{immunoMeta.class}
\alias{immunoMeta.object}
\alias{immunoMeta}
\alias{summary.immunoMeta}
\alias{summary,immunoMeta-method}
\alias{show.immunoMeta}
\alias{show,immunoMeta-method}
\encoding{utf8}

\title{'immunoMeta': a class for storing meta-clustering results}

\description{
The \code{immunoMeta} object contains the clustering results in the 
\emph{immuno}Clust-pipeline obtained by \code{\link{meta.process}}. 
Additionally, it offers methods to structure the meta-clusters and
build up a hierarchical annotation tree. 
}

\usage{

immunoMeta(res,dat,gating)

\S3method{summary}{immunoMeta}(object, \dots)
\S3method{show}{immunoMeta}(object)
}

\arguments{
\item{res}{An \code{immunoClust} object as a result of the meta-clustering.}
\item{dat}{The data on which the meta-clustering was performed.}
\item{gating}{a hierarchial structure annotation of the meta-clusters.}

\item{object}{An object of class \code{immunoMeta} as returned by the 
\code{\link{meta.process}} functions of the 
\emph{immuno}Clust-pipeline.}

\item{\dots}{additinal options for underlying methods.}
}

\value{
An object of class \code{immunoMeta} has the following slots:
\tabular{ll}{
\code{dat.clusters} \tab A \code{dat} list-object of the cell event clusters 
used for meta-clustering.\cr
\code{res.clusters} \tab The 
\code{\link[=immunoClust.object]{immunoClust-object}} of the fitted 
meta-clustering mixture model.\cr
\code{dat.scatter}  \tab A \code{dat} list-object of the scatter parameters for
the cell event clusters used for scatter clustering.\cr
\code{res.scatter}  \tab The 
\code{\link[=immunoClust.object]{immunoClust-object}} of the fitted 
scatter-clustering mixture model.\cr
\code{gating}       \tab A list-object containing the hierarchical 
annotation-tree.\cr
}

The components of the \code{dat} list-objects are:
\tabular{ll}{
\code{P}  \tab The number of parameters for the cell event clusters.\cr
\code{N}  \tab The number of cell-clustering experiments.\cr
\code{K}  \tab The \eqn{N}-dimensional vector with the numbers of cell event 
clusters in each experiment. The total number of clusters is 
\eqn{ sum_{i=1}^N K_i}{sum(K)}.\cr
\code{W}  \tab The \eqn{totK}-dimensional vector with the mixture proportions 
of all clusters.\cr
\code{M}  \tab The \eqn{totK x P}-dimensional matrix of all cluster means.\cr
\code{S}  \tab The \eqn{totK x P x P}-dimensional matrix of all cluster 
covariance matrices.\cr 
\code{expNames}   \tab The \eqn{N}-dimensional character vector with the 
cell-clustering experiment names.\cr
\code{expEvents}  \tab The \eqn{N}-dimensional vector with the numbers of 
events in each cell-clustering experiment.\cr
\code{clsEvents}  \tab The \eqn{totK}-dimensional vector with the number of 
events in each cluster.\cr
\code{desc}     \tab The \eqn{P}-dimensional character vector with the 
parameter description.\cr
}
}


\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}

\seealso{
\code{\link{meta.process}}
}

\examples{
data(dat.meta)
summary(dat.meta)
}
\keyword{print}
