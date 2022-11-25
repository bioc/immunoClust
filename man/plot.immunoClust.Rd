\name{plot.immunoClust}
\docType{methods}
\alias{plot.immunoClust}
\alias{plot,immunoClust-method}
\alias{plot,immunoClust,missing-method}
\alias{plot.immunoClust}

\encoding{utf8}

\title{
Scatterplot of immunoClust Clustering Results
}

\description{  
This method generates scatterplot revealing the cluster assignment.
}

\usage{
\S4method{plot}{immunoClust}(x, data, subset=c(1,2), ellipse=T, 
show.rm=F, include=1:(x@K), main=NULL, 
col=include+1, pch=".", cex=0.6, 
col.rm=1, pch.rm=1, cex.rm=0.6, ecol=col, elty=1, 
npoints=501, add=F, \dots)
}

\arguments{
\item{x}{An object of class \code{\link[=immunoClust.object]{immunoClust}} as 
return by \code{\link{cell.process}}.}
\item{data}{A matrix, data frame of observations, or object of class 
\code{flowFrame}. This is the object of observations on which 
\code{cell.process} was performed or the matrix of cell-cluster centers for 
the \code{meta.process}.}

\item{subset}{A numeric vector of length two indicating which two parameters 
are selected for the scatterplot.  Alternatively, a character vector containing
the names of the two parameters is allowed if \code{x@parameters} is not 
\code{NULL}.}
\item{ellipse}{A logical value indicating whether the cluster 90\% percentil 
boundary is to be drawn or not. }
\item{show.rm}{A logical value indicating whether filtered observations will be
shown or not.}
\item{include}{A numeric vector specifying which clusters will be shown on the 
plot.  By default, all clusters are included.}
\item{main}{Title of the plot.}
\item{col}{Color(s) of the plotting points.  May specify a different color for 
each cluster.}
\item{pch}{Plotting character(s) of the plotting points.  May specify a 
different character for each cluster.} 
\item{cex}{Size of the plotting characters.  May specify a different size for 
each cluster.}
\item{col.rm}{Color of the plotting characters denoting filtered observations.}
\item{pch.rm}{Plotting character used to denote filtered observations.}
\item{cex.rm}{Size of the plotting character used to denote filtered 
observations.}
\item{ecol}{Color(s) of the lines representing the cluster boundaries.  May 
specify a different color for each cluster.}
\item{elty}{Line type(s) drawing the cluster boundaries.  May specify a 
different line type for each cluster.}
\item{npoints}{The number of points used to draw each cluster boundary.}
\item{add}{A logical value.  If \code{TRUE}, add to the current plot.}

\item{\dots}{Further graphical parameters passed to the generic function 
\code{plot}.}
}

\value{
Plots the clustering assignment on an appropriatei plotting device.  
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
dat.res <- dat.exp[[1]]
dat.trans <- trans.ApplyToData(dat.res, dat.fcs)
plot(dat.res, dat=dat.trans,N=1000)
}

\keyword{graphs}
