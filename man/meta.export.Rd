\name{meta.export}
\alias{meta.export}
\alias{meta.numEvents}
\alias{meta.relEvents}
\alias{meta.parMFI}
\alias{meta.numClusters}
\alias{meta.freqTable}
\alias{meta.relEvents2}
\alias{meta.relEvents3}
\alias{meta.majorEvents}
\alias{meta.relParent}

\encoding{utf8}

\title{
immunoClust Meta-clustering Results Export
}

\description{
Thess functions collect the output of the \code{\link{meta.process}} and 
extracts the event numbers, relative frequencies or mean fluorescence 
intensities for each meta-cluster and cell-clustering experiment in a numeric 
table.
}

\usage{
meta.numEvents(meta, out.all=TRUE, out.removed=FALSE, out.unclassified=TRUE)
meta.relEvents(meta, out.all=TRUE, out.removed=FALSE, out.unclassified=TRUE)
meta.relParent(meta, out.all=TRUE, out.removed=FALSE, out.unclassified=TRUE)

meta.parMFI(meta, par, out.all=TRUE)

meta.numClusters(meta, out.all=TRUE)

meta.freqTable(meta)

}

\arguments{
\item{meta}{The list-object returned by the function meta.process.}
\item{par}{An integer index to the specific parameter. }
\item{out.all}{A numeric indicator whether the event numbers of all 
hierarchical gating levels are obtained or only the meta-clusters themselves.}
\item{out.removed}{A numeric indcator whether the number of removed events,
which are not used for clustering are exported.}
\item{out.unclassified}{A numeric indicator whether the event numbers of the 
hierarchical gating levels or all meta-clusters are exported.} 
}

\value{
A numberic matrix with
\describe{
    \item{numEvents}{the number of cell events}
    \item{relEvents}{relative frequencies, i.e. the number of cell events per 
    total meeasured events}
    \item{relParent}{relative frequencies according to parent relationship in
    the annotated hierarchy.}
    \item{parMFI}{mean fluorecence intensities in one parameter, i.e. the 
    meta-cluster centers in asinh-tranformed scale}
    \item{numClusters}{the number of cell clusters}
    \item{freqTable}{relative frequencies with respect to all gating 
    hierarchie levels}
} 
in each meta-cluster (and gating hierarchy level) for each cell-clustering 
experiment.
}

\references{
Sörensen, T., Baumgart, S., Durek, P., Grützkau, A. and Häupl, T. 
immunoClust - an automated analysis pipeline for the identification of 
immunophenotypic signatures in high-dimensional cytometric datasets. 
\emph{Cytometry A} (submitted). 
}

\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}

\seealso{
\code{\link{meta.process}}
}

\examples{
data(dat.exp)
meta <- meta.process(dat.exp)
tbl <- meta.numEvents(meta)
}
