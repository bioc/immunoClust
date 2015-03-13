\name{meta.exprs}
\alias{meta.exprs}

\encoding{utf8}

\title{
Collecting Data of an immunoClust vector
}
\description{
The function takes a vector of \code{immunoClust}-object obtained
by the \code{cell.process} function and extracts ths information
into a \code{list} object.
}
\usage{
meta.exprs(exp, sub=c())
}

\arguments{
\item{exp}{The vector of \code{immunoClust} object with the cell clustering 
results.}
\item{sub}{A integer array indicating the parameter subset to be collected.}
}

\value{
A \code{list} object with the following slots:

\item{P}{The number of observed parameters for the cell event clusters.}
\item{N}{The number of cell-clustering samples.}
\item{K}{The \eqn{N}-dimensional vector with the numbers of cell event 
clusters in each experiment. The total number of clusters is  
\eqn{totK = sum_{i=1}^N K_i}{sum(K)}.}
\item{W}{The \eqn{totK}-dimensional vector with weights of all clusters.}
\item{M}{The \eqn{totK x P}-dimensional matrix of all cluster means.}
\item{S}{The \eqn{totK x P x P}-dimensional matrix of all cluster covariance 
matrices.}
\item{expNames}{The \eqn{N}-dimensional vector with the experiment names of
the cell clustering samples.}
\item{expEvents}{The \eqn{N}-dimensional vector for the total number of events
of the cell clustering samples.}
\item{clsEvents}{The \eqn{totK}-dimensional vector for the event number of all
clusters.}
\item{desc}{The \eqn{P}-dimensional vector for the parameter description.} 
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
\code{\link[=immunoClust.object]{immunoClust}}.
}
\examples{
data(dat.exp)
d <- meta.exprs(dat.exp, sub=c(1,2))
}
