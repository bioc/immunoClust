\name{trans.ApplyToData}
\alias{trans.ApplyToData}

\encoding{utf8}

\title{
immunoClust asinh-Transformation
}

\description{
Applies the transformation information of the \code{immunoClust} object to the 
raw observed FC dataset.
}

\usage{
trans.ApplyToData(x, data, max.decade=attr(x,"trans.decade"), 
    lin.scale=attr(x,"trans.scale") )
}

\arguments{
\item{x}{The \code{immunoClust} object containing the estimators for the 
transformation \code{trans.a} and \code{trans.b}.}
\item{data}{The numeric matrix, data frame of observations, or object of class 
\code{flowFrame}.}
\item{max.decade}{A numeric scale for the maximum transformed observation value;
if missing or below 0, no scaling of the transformed values is apllied, which 
is the default in \emph{immuno}Clust.}
\item{lin.scale}{A numeric scaling factor for the linear, i.e. not transformed,
parameters; if missing no scaling, i.e. \eqn{lin.scale=1}, is applied, which is
the default in \emph{immuno}Clust.}
}

\details{
In \emph{immuno}Clust an \eqn{asinh}-transformation 
\eqn{h(y)=asinh(a \cdot y + b)}{h(y)=asinh(a*y + b)} is applied to the 
fluorescence parameter in the observed data. The scatter parameter are assumed 
to be linear. 
}

\value{
A matrix or \code{flowFrame} with replaced transformed oberservation values.
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
\code{\link[=immunoClust.object]{immunoClust}}, 
\code{\link{trans.FitToData}}, \code{\link{cell.process}}
}

\examples{
data(dat.fcs)
data(dat.exp)
dat.trans <- trans.ApplyToData(dat.exp[[1]], dat.fcs)
\dontrun{
plot(dat.exp[[1]], data=dat.trans)
}
}
