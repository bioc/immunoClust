\name{cell.EM}
\alias{cell.EM}
\alias{cell.Estimation}

\encoding{utf8}

\title{
immunoClust EMt-iteration on Cell-events given initial Model Parameters
}

\description{
Performs EMt-iteration on cell event observations giving initial model 
parameters and returns the fitted clusters information in an object of class 
\code{\link[=immunoClust.object]{immunoClust}}.
}

\usage{
cell.EM(data, parameters=NULL, expName="immunoClust Experiment", 
    history=NULL, state=NULL,
    K, w, m, s, B=50, tol=1e-5, bias=0.5, modelName="mvt")

cell.Estimation(data, parameters=NULL, expName="immunoClust Experiment",  
    history=NULL, state=NULL,
    K, w, m, s, scale_Z=TRUE, modelName="mvt")
}

\arguments{
\item{data}{A numeric matrix, data frame of observations, or object of class 
flowFrame.} 
\item{parameters}{ A character vector specifying the parameters (columns) to be
included in clustering. When it is left unspecified, all the parameters will be
used.}
\item{expName}{The name of the clustering experiment.}
\item{history}{experimental; unused so far.}
\item{state}{experimental: unused so far.}
\item{K}{The number of clusters. }
\item{w}{The \eqn{K}-dimensional vector of the mixture proportions.}
\item{m}{The \eqn{K x P}-dimensional matrix of the \eqn{K} estimated cluster 
means. }
\item{s}{The \eqn{K x P x P}-dimensional matrix of the \eqn{K} estimated cluster
covariance matrices.}
\item{B}{The maximum number of EMt-iterations.}
\item{tol}{The tolerance used to assess the convergence of the EMt-algorithms.}
\item{bias}{The ICL-bias used in the EMt-algorithm.}
\item{scale_Z}{Scale the returned a-posteriori probabilities to one for each observed event.}
\item{modelName}{Used mixture model; either \code{"mvt"} or \code{"mvn"} for a 
\eqn{t}- or Gaussian mixture model respectively.}
}

\details{
Whereas \code{cell.EM} performs a complete EMt-iteration, \code{cell.Estimate} 
only calculates the a-posteriori probabilities and the Maximum-A-Posteriori 
estimators of cluster membership for all events.
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
\code{\link{cell.ME}}, \code{\link{cell.FitModel}}
}
\examples{
data(dat.fcs)
data(dat.exp)
## cell.clustering result for dat.fcs
r <- dat.exp[[1]]
summary(r)
## apply model parameter to all (unfiltered) events
dat.trans <- trans.ApplyToData(r, dat.fcs)
r2 <- cell.EM(dat.trans, parameters(r), K=ncls(r), 
            w=weights(r),m=mu(r),s=sigma(r))
summary(r2)
}
\keyword{ cluster }
