\name{immunoClust-object}
\alias{immunoClust.object}
\alias{immunoClust-class}
\alias{summary}
\alias{summary,immunoClust-method}
\alias{summary.immunoClust}
\alias{show}
\alias{show,immunoClust-method}
\alias{show.immunoClust}

\encoding{utf8}

\title{immunoClust-Object}

\description{
The \code{immunoClust} object contains the clustering results in the 
\emph{immuno}Clust-pipeline as obtained by \code{\link{cell.process}} or 
\code{\link{meta.process}}.
}

\usage{
\S4method{summary}{immunoClust}(object)
\S4method{show}{immunoClust}(object)
}

\arguments{
\item{object}{An object of class \code{immunoClust} as returned by the 
\code{\link{cell.process}} or \code{\link{meta.process}} functions of the 
\emph{immuno}Clust-pipeline.}
}

\value{
An object of class \code{immunoClust} has the following slots:

\tabular{ll}{
\code{expName} \tab The name of the clustering experiment.\cr
\code{fcsName} \tab The path of the clustered FCS-file.\cr
\code{parameters} \tab The parameters used for clustering.\cr
\code{removed.below} \tab Number of observations removed from below.\cr
\code{removed.above} \tab Number of observations removed from above.\cr
\code{trans.a} \tab The \eqn{P}-dimensional vector of the scaling factors for 
the asinh-transformation of each used parameter. A scaling factor of 0 indicates
that a parameter is not transformed.\cr
\code{trans.b} \tab The \eqn{P}-dimensional vector of the translations for the 
asinh-transformation of each used parameter.\cr
\code{trans.decade} \tab experimental; should be -1.\cr
\code{trans.scale} \tab experimental; should be 1.0.\cr
\code{K} \tab The number of clusters.\cr
\code{N} \tab The number of observations.\cr
\code{P} \tab The number of used parameters.\cr
\code{w} \tab The \eqn{K}-dimensional vector of the mixture proportions.\cr
\code{mu} \tab The \eqn{K x P}-dimensional matrix of the \eqn{K} estimated 
cluster means.\cr
\code{sigma} \tab The \eqn{K x P x P}-dimensional matrix of the \eqn{K} 
estimated cluster covariance matrices.\cr
\code{z} \tab The \eqn{K x N}-dimensional matrix containing the a-posteriori 
probabilities of cluster membership.\cr
\code{label} \tab The \eqn{N}-dimensional vector containing the maximum a 
posteriori estimator for cluster membership.\cr
\code{logLike} \tab A vector of length 3 containing the BIC, ICL and the 
classification likelihood without penalty of the fitted model.\cr
\code{BIC} \tab The Bayesian Information Criterion for the fitted mixture 
model.\cr
\code{ICL} \tab The Integrate Classification Likelihood for the fitted model.\cr
\code{history} \tab experimental; unused so far.\cr
\code{state} \tab experimental; unused so far.\cr
}
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
\code{\link{cell.process}}, \code{\link{meta.process}}
}

\examples{
data(dat.exp)
summary(dat.exp[[1]])
}
\keyword{print}
