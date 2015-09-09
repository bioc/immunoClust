\name{dat.meta}
\alias{dat.meta}
\docType{data}

\title{
immunoClust Meta-clustering Results Sample
}
\description{
The Meta-clustering result of the \code{\link{dat.exp}} data set.
}

\usage{data("dat.meta")}

\value{
A list-object containing the meta-clusering result. A detailed 
description is documented in the value section for the
\code{\link{meta.process}} function.
}

\details{
The Meta-clustering was performed with an ICL-bias of 0.4.
}

\source{
http://flowrepository.org/id/FR-FCM-ZZWB
}

\examples{
data(dat.meta)

## extract event counts in the 5 samples for all meta clusters
res <- meta.numEvents(dat.meta)

}
\keyword{datasets}
