\name{generics.immunoclust}
\docType{class}
\alias{nsam}
\alias{sam_ncls}
\alias{sam_clsWeights}
\alias{sam_clsEvents}
\alias{sam_clsMu}
\alias{sam_clsSigma}
\alias{nobs}
\alias{npar}
\alias{ncls}
\alias{weights}
\alias{mu}
\alias{sigma}
\alias{label}
\alias{posterior}
\alias{subset}
\alias{parameters}
\alias{transformParams}

\encoding{utf8}

\title{
Generic function definitions in immunoClust
}

\description{
Collection of generic function definitions used in immunoClust 
either for an \code{immunoClust} or an \code{immunoMeta} object.
}

\usage{
nsam(object, \dots)

sam_ncls(object, \dots)

sam_clsWeights(object, \dots)

sam_clsEvents(object, \dots)

sam_clsMu(object, \dots)

sam_clsSigma(object, \dots)

nobs(object, \dots)

npar(object, \dots)

ncls(object, \dots)

weights(object, \dots)

mu(object, \dots)

sigma(object, \dots)

label(object, \dots)

posterior(object, \dots)

subset(x, \dots)

parameters(object, \dots)

transformParams(object, \dots)
}

\arguments{
\item{object, x}{an object to apply the function.}
\item{\dots}{addionional options to be passed tomethods}
}

\value{
The appropriate value for the specific cal (see dection Details).
}

\section{Details}{
\describe{
    \item{nsam}{returns the number of cell-event \code{immunoClust}-objects 
    co-clustered in the \code{immunoMeta}-object.}

    \item{sam_clsWeights}{returns the cluster weights of all samples 
    cell-clusters.}

    \item{sam_clsEVents}{returns the cluster event numbers of all samples 
    cell-clusters.}

    \item{sam_clsMu}{returns the cluster means of all samples cell-clusters.}

    \item{sam_clsSigma}{returns the cluster co-variance matrices of all samples 
    cell-clusters.}

    \item{nobs}{already generic in \code{stats}. Here, returns the number of 
    clustered objects either cell-events or cell-clusters in cell event or 
    meta clustering.}

    \item{npar}{returns the number of parameters used for clustering.}

    \item{ncls}{returns the number of cluster, either cell-event cluster or 
    meta-cluster.}

    \item{weights}{already generic in \code{stats}. Here, returns the weights 
    of the mixture models for the cell-event of meta-clustering.}

    \item{mu}{returns the cluster means.}

    \item{sigma}{already generic in \code{stats}. Here, returns the co-variance
    matrices of the clusters.}

    \item{label}{returns the cluster label, i.e. the assignment of the 
    clustered objects to the clusters.}
    
    \item{posterior}{returns the posterior probabilities of cluster membership
    for the clustered objects.}

    \item{events}{returns the number of cell-events for the clusters.}

    \item{subset}{alreay generic in \code{stats}. Here, returns an object with
    mixture model on a subset of parameters.}

    \item{parameters}{already generic in \code{flowCore}. Here, lists the
    parameters used for clustering.}

    \item{transformParam}{return an object with transformed mixture model
    parameters.}
}
}



\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}

\seealso{
\code{\link{immunoClust}}, \code{\link{immunoMeta}}
}


\keyword{ generics }
