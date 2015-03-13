\name{trans.FitToData}
\alias{trans.FitToData}

\encoding{utf8}

\title{
immunoClust asinh-Transformation Optimization
}

\description{
Performs variance stabilization transformation estimation on the fluorescense 
parameters of the observed cell events. It is integrated in the interative 
cell event clustering approach of \emph{immuno}Clust when transformation 
estimation should be applied. 
}

\usage{
trans.FitToData(x, data, B=10, tol=1e-5, certainty=0.3, proc="vsHtransAw")
}

\arguments{
\item{x}{The \code{immunoClust} object of the fitted mixture model and initial 
estimators for the transformation.}
\item{data}{The numeric matrix, data frame of observations, or object of class 
flowFrame.}
\item{B}{The maximum number of BFG2 minimizer iterations.}
\item{tol}{The tolerance used to assess the convergence for the BFG2 minimizer.}
\item{certainty}{Minimum probability for cluster membership of an observation 
to be taken into account.}
\item{proc}{An experimental switch for alternative procedures; should be 
"vsHtransAw".}
}

\details{
In \emph{immuno}Clust an \eqn{asinh}-transformation 
\eqn{h(y)=asinh(a \cdot y + b)}{h(y)=asinh(a * y + b)} is applied for all 
fluorescence parameter in the observed data. 

The transformation optimization \code{trans.FitToData} requires a fitted model 
of cluster information together with suitable initial transformation estimation
in an \code{\link[=immunoClust.object]{immunoClust}} object. It fits the 
transformation based on the initial scaling values \code{trans.a} and 
translation values \code{trans.b} to the observed data. 
It returns the optimized transformation parameter in a 
\eqn{2 \times P}{2 x P}-dimensional matrix, first row for the scaling and 
second row for the translation values. 
A scaling value of \eqn{a=0} on input and output indicates, that a parameter 
should not be transformed. 

The presented transformation optimization ("vsHtransAw") fits only the scaling 
value. An alternative procedure ("vsHtrans_w") fits both, the scaling and the 
translation value, but turns out to be less robust.
}

\value{
Optimized transformation scaling and translation values in a 
\eqn{2 \times P}{2 x P}-dimensional matrix, first row for the scaling and 
second row for the translation values.
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
\code{\link{trans.ApplyToData}}, \code{\link{cell.process}}
}

\examples{
data(dat.fcs)
data(dat.exp)
## in dat.exp the z-matrices of the immunoClust-object are removed
## so we have to re-calculate it first ...
dat.trans <- trans.ApplyToData(dat.exp[[1]], dat.fcs)
res <- cell.Classify(dat.exp[[1]], dat.trans)
## ... now the transformation parameter can be optimzed
trans.FitToData(res, dat.fcs)
}
