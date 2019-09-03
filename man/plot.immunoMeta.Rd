\name{plot.immunoMeta}
\docType{methods}
\alias{plot.immunoMeta}
\alias{plot,immunoMeta-method}
\alias{plot,immunoMeta,missing-method}


\encoding{utf8}

\title{
Scatterplot of immunoMeta Structured Clustering Results
}

\description{  
This method generates scatterplot revealing the cluster assignment.
}

\usage{

\method{plot}{immunoMeta}(x, pos=c(), main="", plot.childs=TRUE,
plot.unclassified=FALSE, plot.subset=c(), inc.childs=c(), plot.ellipse=TRUE,
plot.all=FALSE, \dots)

}

\arguments{
\item{x}{An object of class \code{\link[=immunoMeta.object]{immunoMeta}} as 
return by \code{\link{meta.process}}.}

\item{pos}{gives the position in the immunoMeta-hierarchy to plot
(\code{default=c()} plots the root level). 
\code{pos} is an array of indices, which addresses the level of interest.
Each level in the immunoMeta-hierarchy has an array of sub-levels and
\code{pos} is the sequences of indices into these sub-levels.} 

\item{main}{additional title which is concatenated with the position and
description path of the plotted level.}

\item{plot.subset}{an array of indices for the parameter selection to be
plotted.}

\item{plot.unclassified}{if set, the unclassified  clusters,i.e clusters not
assigned into a sub-level, are plotted rather than the classified clusters.}

\item{plot.childs}{colours the clusters by the sub-level rather than the
clusters themselves. By default colours are assigned by sub-level index
repeated in red, green,blue,cyan,magenta,yellow,gray,black}

\item{inc.childs}{optionally, to restrict to a particular selection of
sub-levels to plot.}

\item{plot.ellipse}{surrounds the cell-cluster center by an ellipse reflecting
the meta-cluster deviation}

\item{plot.all}{plots all sub-levels. Usefull for a full annotation
documentation with a pdf file.}

\item{\dots}{Further graphical parameters passed to the generic function
\code{plot}.}
}

\value{
Plots the clustering assignment on an appropriated plotting device.  
}

\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}


\seealso{
\code{\link[=immunoMeta.class]{immunoMeta-object}}
}

\examples{
data(dat.meta)
plot(dat.meta)
}

\keyword{graphs}
