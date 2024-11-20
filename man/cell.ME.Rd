\name{cell.ME}
\alias{cell.MEstep}
\alias{cell.Mstep}

\encoding{utf8}

\title{
immunoClust EM-iteration on Cell-events given initial Cluster Membership 
Assignment
}

\description{
Performs an EM-iteration on cell event observations given an initial cluster 
membership for the cell events and returns the fitted cluster information in an
object of class \code{\link[=immunoClust.object]{immunoClust}}.
}

\usage{

cell.MEstep(data, label, parameters=NULL, 
    expName="immunoClust Experiment", 
    B=1, tol=1e-5, modelName="mvt")
        
cell.Mstep(data, label, parameters=NULL,
    expName="immunoClust Mstep", modelName="mvt")
}

\arguments{
\item{data}{A numeric matrix, data frame of observations, or object of class 
flowFrame.} 
\item{parameters}{ A character vector specifying the parameters (columns) to be
included in clustering. When it is left unspecified, all the parameters will be
used.}
\item{expName}{The name of the clustering experiment.}
\item{label}{ The \eqn{N}-dimensional vector containing the initial cluster 
membership. A label-number of 0 for an event indicates that this event is not 
initially assigned to a cluster.}
\item{B}{The maximum number of EMt-iterations.}
\item{tol}{The tolerance used to assess the convergence of the EMt-algorithms.}
\item{modelName}{Used mixture model; either \code{"mvt"} or \code{"mvn"} for a 
\eqn{t}- or Gaussian mixture model respectively.}
}

\details{
\code{cell.ME} and \code{cell.MEstep} do the same call. In \code{cell.MEstep}
the calling options are a bit better structured and \code{cell.ME}
becomes deprecated in future.
}
\value{
The fitted clusters information in an object of class 
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
\code{\link{cell.EMt}}
}

\examples{
data(dat.fcs)
data(dat.exp)
## cell.clustering result for dat.fcs
r1 <- dat.exp[[1]]
summary(r1)
## apply model parameter to all (unfiltered) events
dat.trans <- trans.ApplyToData(r1, dat.fcs)
r2 <- cell.MEstep(dat.trans, label(r1), parameters=parameters(r1) )
summary(r2)
}
\keyword{ cluster }
