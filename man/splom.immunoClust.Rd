\name{splom.immunoClust}
\docType{methods}
\alias{splom}
\alias{splom.immunoClust}
\alias{splom,immunoClust,missing-method}
\alias{splom,immunoClust,flowFrame-method}
\alias{splom,immunoClust,matrix-method}
\alias{datSplom}

\encoding{utf8}

\title{
Scatterplot Matrix of immunoClust Clustering Results
}

\description{
This method generates scatterplot matrix revealing the cluster assignment.
}

\usage{
\S4method{splom}{immunoClust,missing}(x, data, include=seq_len(x@K), \dots)

\S4method{splom}{immunoClust,flowFrame}(x, data, include=seq_len(x@K), 
subset=seq_len(length(attributes(x)$param)), N=NULL,label=NULL, desc=NULL, 
add.param=c(), \dots)

\S4method{splom}{immunoClust,matrix}(x, data, include=seq_len(x@K), 
subset=seq_len(length(attributes(x)$param)), N=NULL, label=NULL, 
desc=NULL, \dots)

datSplom(label, data, subset=seq_len(ncol(data)), 
include=seq_len(nrow(data)), \dots) 
}

\arguments{
\item{x}{An object of class \code{\link[=immunoClust.object]{immunoClust}} as 
return by \code{\link{cell.process}} or \code{\link{meta.process}}.}
\item{data}{Missing, a matrix, or object of class \code{flowFrame}. This is the
object of observations on which \code{cell.process} was performed.}
\item{include}{A numeric vector specifying which clusters will be shown on the 
plot.  By default, all clusters are included.}
\item{subset}{A numeric vector indicating which parameters are selected for the
scatterplot matrix.} 
\item{N}{An integer for the maximum number of observations to be plotted. By 
default all observations are plotted.} 
\item{label}{A integer vector for the cluster mebership of the observations. By
default this is \code{x@label}.}
\item{desc}{A character vector for the parameter description. }
\item{add.param}{ A list of additional parameters to plot, which are not used 
for clustering.}
\item{\dots}{Further graphical parameters passed to the generic function 
\code{splom}.}
}

\value{
An object of class \code{trellis} as returned by the generic \code{splom} 
function of the \code{lattice}-package. The \code{print} method (called by 
default) will plot it on an appropriate plotting device.
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
\code{\link[=immunoClust.object]{immunoClust-object}}
}
\examples{
data(dat.fcs)
data(dat.exp)
# cell clustering results of dat.fcs
dat.res <- dat.exp[[1]]
dat.trans <- trans.ApplyToData(dat.res, dat.fcs)
splom(dat.res, data=dat.trans, N=1000)
}

\keyword{graphs}
