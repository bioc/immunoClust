\name{cell.hclust}
\alias{cell.hclust}

\encoding{utf8}

\title{
Hierarchical Model Based Clustering of Cell-events in the immunoClust-pipeline
}

\description{
Performs model based agglomerative clustering on cell event observations with 
weights. It is used in the interative cell event clustering approach of 
\emph{immuno}Clust to obtain an initial cluster membership for the 
EM(t)-iteration. 
}

\usage{
cell.hclust(data, weights=NULL)
}

\arguments{
\item{data}{The numeric \eqn{N \times P}{N x P}-dimensional data matrix to 
cluster. Each row contains a \eqn{P}-dimensional overservation vector.}
\item{weights}{The \eqn{N}-dimensional vector of optional weights to be applied
for the overservations.}
}
\details{
This function is used internally in \code{\link{cell.TestSubCluster}} procedure
of \pkg{immunoClust}.
}

\value{
A numeric \eqn{(N-1) \times 2}{(N-1) x 2}-dimensional matrix which gives the 
minimum index for observations in each of the two clusters merged at the 
\emph{i}th step in each row.
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
\code{\link{cell.TestSubCluster}}, \code{\link{cell.process}}
}

\examples{
data(dat.fcs)
inc <- sample(1:nrow(dat.fcs), 50)
result <- cell.hclust(exprs(dat.fcs)[inc,])
}

\keyword{ cluster }
