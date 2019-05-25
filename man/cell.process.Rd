\name{cell.process}
\alias{cell.process}
\alias{cell.MajorIterationLoop}
\alias{cell.MajorIterationTrans}
\alias{cell.InitialModel}
\alias{cell.classifyAll}

\encoding{utf8}

\title{
Clustering of Cell-events in the immunoClust-pipeline
}

\description{
This function performs iterative model based clustering on cell-event data. It 
takes the observed cell-event data as major input and returns an object of class
\code{immunoClust}, which contains the fitted mixture model parameter and 
cluster membership information. The additional arguments control the routines 
for data preprocessing, major loop and EMt-iteration, the model refinement 
routine and transformation estimation.
}

\usage{
cell.process(fcs, parameters=NULL, 
    apply.compensation=FALSE, classify.all=FALSE, 
    N=NULL, min.count=10, max.count=10, min=NULL, max=NULL,  
    I.buildup=6, I.final=4, I.trans=I.buildup, 
    modelName="mvt", tol=1e-5, bias=0.3,
    sub.tol= 1e-4, sub.bias=bias, sub.thres=bias, sub.samples=1500, 
    sub.extract=0.8, sub.weights=1, sub.standardize=TRUE,
    trans.estimate=TRUE, trans.minclust=10, 
    trans.a=0.01, trans.b=0.0, trans.parameters=NULL)

cell.MajorIterationLoop(dat, x=NULL, parameters=NULL, 
    I.buildup=6, I.final=4, 
    modelName="mvt", tol=1e-5, bias=0.3,
    sub.bias=bias, sub.thres=0.0, sub.tol=1e-4, sub.samples=1500, 
    sub.extract=0.8, sub.weights=1, sub.EM="MEt", sub.standardize=TRUE, seed=1)

cell.MajorIterationTrans(fcs, x=NULL, parameters=NULL, 
    I.buildup=6, I.final=4, I.trans=I.buildup, 
    modelName="mvt", tol=1e-5, bias=0.3,
    sub.bias=bias, sub.thres=0.0, sub.tol=1e-4, sub.samples=1500, 
    sub.extract=0.8, sub.weights=1, sub.EM="MEt", sub.standardize=TRUE, seed=1, 
    trans.minclust=5, trans.a=0.01, trans.decade=-1, trans.scale=1.0, 
    trans.proc="vsHtransAw")

cell.InitialModel(dat, parameters=NULL, trans.a = 0.01, trans.b = 0.0, 
    trans.decade=-1, trans.scale=1.0)

cell.classifyAll(fcs, x, apply.compensation=FALSE)                         
}

\arguments{
\item{fcs}{An object of class flowFrame. Rows correspond to observations and 
columns correspond to measured parameters.}
\item{dat}{A numeric matrix, data frame of observations, or object of class 
flowFrame. Rows correspond to observations and columns correspond to measured 
parameters.}
\item{x}{An object of class \code{immunoClust}. Used as inital model int the 
major iteration loop. When left unspecified the simplest model containing 1 
cluster is used as initial model.}

\bold{ Arguments for data pre and post processing: }
\item{parameters}{A character vector specifying the parameters (columns) to be 
included in clustering. When it is left unspecified, all the parameters will be
used.}
\item{apply.compensation}{A numeric indicator whether the compensation matrix 
in the flowFrame should be applied.}
\item{classify.all}{A numeric indicator whether the removed over- and 
underexposed observations should also be classified after the clustering 
process.}
\item{N}{Maximum number of observations used for clustering. When unspecified 
or higher than the number of observations (i.e. rows) in dat, all observations 
are used for clustering, otherwise only the first \code{N} observations.}
\item{min.count}{An integer specifying the threshold count for filtering data 
points from below. The default is 10, meaning that if 10 or more data points 
are smaller than or equal to \code{min}, they will be excluded from the 
analysis. If \code{min} is \code{NULL}, then the minimum value of each 
parameter will be used. To suppress filtering, it is set to -1.}
\item{max.count}{An integer specifying the threshold count for filtering 
data points from above. Interpretation is similar to that of \code{min.count}.}
\item{min}{The lower limit set for data filtering. Note that it is a vector of 
length equal to the number of parameters (columns), implying that a different 
value can be set for each parameter.}
\item{max}{The upper limit set for data filtering. Interpretation is similar to
that of \code{min}.}

\bold{Arguments for the major loop and EMt-iteration:}
\item{I.buildup}{The number of major iterations, where the number of used 
observations is doubled successively.}
\item{I.final}{The number of major iterations with all observations.}
\item{I.trans}{The number of iterations where transformation estimation is 
applied.}
\item{modelName}{Used mixture model; either \code{"mvt"} for a \eqn{t}-mixture 
model or \code{"mvn"} for a Gaussian Mixture model. With \code{"mvt2"} an 
implementation variant for \code{"mvt"} is given, which is more reliable for 
samples with cutted values at the lower or upper edges of the parameter space 
(e.g. for CyTOF all values below a detection limit are set to zero which leads 
to wrong co-variance estimators and poor clustering results).}
\item{tol}{The tolerance used to assess the convergence of the major 
EM(t)-algorithms of all observations.}
\item{bias}{The ICL-bias used in the major EMt-algorithms of all observations.}

\bold{Arguments for model refinement (sub-clustering):}
\item{sub.tol}{The tolerance used to assess the convergence of the 
EM-algorithms in the sub-clustering.}
\item{sub.bias}{The ICL-bias used in the sub-clustering EMt-algorithms, in 
general the same as the ICL-bias.}
\item{sub.thres}{Defines the threshold, below which an ICL-increase is 
meaningless. The threshold is given as the multiple (or fraction) of the costs 
for a single cluster.}
\item{sub.samples}{The number of samples used for initial hierarchical 
clustering.}
\item{sub.extract}{The threshold used for cluster data extraction.}
\item{sub.weights}{Power of weights applied to hierarchical clustering, where 
the used weights are the probabilities of cluster membership.}
\item{sub.EM}{Used EM-algorithm; either \code{"MEt"} for EMt-iteration or 
\code{"ME"} for EM-iteration without test step.}
\item{sub.standardize}{A numeric indicating whether the samples for 
hierarchical clustering are standardized (mean=0, SD=1).}
\item{seed}{The seed integer for the random number generator.}

\bold{Arguments for transformation optimization:}
\item{trans.estimate}{A numeric indicator whether transformation estimation 
should be applied.}
\item{trans.minclust}{The minimum number of clusters required to start 
transformation estimation.}
\item{trans.a}{A numeric vector, giving the (initial) scaling \eqn{a} for the 
asinh-transformation \eqn{h(y) = asin(a \cdot y + b)}. A scaling factor of 
\eqn{a=0} indicates that a parameter is not transformed.}
\item{trans.b}{A numeric vector, giving the (initial) translation \eqn{b} for 
the asinh-transformation.}
\item{trans.parameters}{A character vector, specifying the parameters (columns)
to be applied for transformation. When it is left unspecified, the parameters 
to be transformed are obtained by the \code{PxDISPLAY} information of the 
flowFrame description parameters. All parameters with LOG display values are 
transformed.}
\item{trans.decade}{A numeric scale value for the theorectical maximum of 
transformed observation value. If below 0, no scaling of the trasnformed values
is applied, which is the default in the \emph{immuno}Clust-pipeline.}
\item{trans.scale}{A numeric scaling factor for the linear (i.e. not 
transformed) parameters. By default the linear parameters (normally the scatter
parameters) are not scaled.}
\item{trans.proc}{An experimental switch for alternative procedures; should 
be "vsHtransAw".}
}

\details{
The \code{cell.process} function does data preprocessing and calls the major 
iteration loop either with or without integrated transformation optimization.
When transformation optimization is applied the transformation parameters give 
the initial transformation otherwise they define the fixed
transformation.

The major iteration loop with included transformation optimization relies on 
\code{flowFrames} structure from the \code{flowCore}-package for the storage of
the observed data.

The \code{cell.InitialModel} builds up an initial \code{immunoClust-object} 
with one cluster and the given transformation parameters.

The \code{cell.classifyAll} calculates the cluster membership for the removed 
cell events. The assigment of the cluster membership is critical for over- and 
underexposed obsevervations and the interpretaion is problematic.
}

\value{
The fitted model information in an object of class 
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
\note{
a) The data preprocessing arguments (\code{min.count}, \code{max.count}, 
\code{min} and \code{max}) for removing over- and underexposed observations are
adopted from \code{flowCust-package} with the same meaning.

b) The \code{sub.thres} value is given in here in relation to the single 
cluster costs 
\eqn{\frac{1}{2} \cdot P \cdot (P+1) \cdot log(N)}{1/2 x P x (P+1) x log(N)}. 
An absolute increase of the log-likelihood above is reported as 
reasonable from the literature. From our experience a higher value is required 
for this increase in FC data. For the ICL-bias and the sub.thres identical 
values were chosen. For the CyTOF dataset this value had been adjusted to 0.05 
since the absolute increase of the log-likelihood became to high due to the 
high number of parameters.

c) The \code{sub.extract} value controls the smooth data extraction for a 
cluster. A higher value includes more events for a cluster in the 
sub-clustering routine. 

d) The default value of \code{trans.a=0.01} for the initial transformation is 
optimized for Fluorescence Cytometry. For CyTOF data the initial scaling value 
was \code{trans.a=1.0}. 
}

\seealso{
\code{\link[=immunoClust.object]{immunoClust-object}}, 
\code{\link[=plot.immunoClust]{plot}},
\code{\link[=splom.immunoClust]{splom}},
\code{\link{cell.FitModel}},
\code{\link{cell.SubClustering}},
\code{\link{trans.FitToData}}
}

\examples{
data(dat.fcs)
res <- cell.process(dat.fcs)
summary(res)
}

\keyword{ cluster }
