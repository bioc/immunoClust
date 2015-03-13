\name{dat.exp}
\alias{dat.exp}
\docType{data}

\title{
immunoClust Meta-clustering Sample
}
\description{
A \code{vector} of \code{immunoClust-objects} with \code{cell.process} 
clustering results of five samples. 
}

\usage{data("dat.exp")}

\value{
A vector of 5 \code{\link[=immunoClust.object]{immnoClust-objects}} for the 
cell clustering results of 5 FC samples.
\describe{
    \item{\code{[[1]]}}{CD19 MACS-depleted cells}
    \item{\code{[[2]]}}{CD15 MACS-depleted cells}
    \item{\code{[[3]]}}{CD14 MACS-depleted cells}
    \item{\code{[[4]]}}{CD4 MACS-depleted cells}
    \item{\code{[[4]]}}{CD3 MACS-depleted cells}
}
}

\details{
Cell-event clustering was performed on reduced (10.000 events) sample data of 
the dataset of \emph{immuno}Clust, MACS-depleted populations datasets 2010. 
URL http://flowrepository.org/id/FR-FCM-ZZWB.
}

\source{
http://flowrepository.org/id/FR-FCM-ZZWB
}

\examples{
data(dat.exp)

## process meta clustering
meta <- meta.process(dat.exp, meta.bias=0.6)

## extract event counts in the 5 samples for all meta clusters
res <- meta.numEvents(meta)

}
\keyword{datasets}
