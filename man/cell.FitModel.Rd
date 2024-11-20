\name{cell.FitModel}
\alias{cell.FitModel}
\alias{cell.Classify}

\encoding{utf8}

\title{
immunoClust EMt-iteration on Cell-events given initial Model Parameters
}
\description{
The function fits initial model parameters to specific observed cell event data.
The function returns the cluster information of the fitted model in an object 
of class \code{\link[=immunoClust.object]{immunoClust}}.
}
\usage{
cell.FitModel(x, data, B=50, tol=1e-5, bias=0.5, modelName="mvt" )

cell.Classify(x, data, modelName="mvt" )
}

\arguments{
\item{x}{An immunoClust object with the initial model parameter 
(\eqn{parameters}, \eqn{K}, \eqn{w}, \eqn{mu}, \eqn{sigma}).}
\item{data}{A numeric matrix, data frame of observations, or object of class 
flowFrame.} 
\item{B}{The maximum number of EMt-iterations.}
\item{tol}{The tolerance used to assess the convergence of the EMt-algorithms.}
\item{bias}{The ICL-bias used in the EMt-algorithm.}
\item{modelName}{Used mixture model; either \code{"mvt"} or \code{"mvn"} for a 
\eqn{t}- or Gaussian mixture model respectively.}
}
\details{
These functions are wrapper of the functions \code{cell.EM} and 
\code{cell.Estimation}, when model cluster parameters are combined in an object
of class \code{immunoClust} and are used in the iterative cell event clustering
process \code{\link{cell.process}} of \emph{immuno}Clust and are not intended to
be called directly.
}

\value{
The fitted model cluster information in an object of class 
\code{\link[=immunoClust.object]{immunoClust}}.
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
\code{\link{cell.process}}, \code{\link{cell.EMt}}, \code{\link{cell.Estep}}
}

\examples{
data(dat.fcs)
data(dat.exp)
r1 <- dat.exp[[1]]
dat.trans <- trans.ApplyToData(r1, dat.fcs)
r2 <- cell.FitModel(r1, dat.trans)
}
