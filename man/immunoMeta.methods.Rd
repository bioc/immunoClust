\name{methods.immunoMeta}
\docType{class}
\alias{nsam.immunoMeta}
\alias{nsam,immunoMeta-method}
\alias{sam_ncls.immunoMeta}
\alias{sam_ncls,immunoMeta-method}
\alias{sam_clsWeights.immunoMeta}
\alias{sam_clsWeights,immunoMeta-method}
\alias{sam_clsMu.immunoMeta}
\alias{sam_clsMu,immunoMeta-method}
\alias{sam_clsSigma.immunoMeta}
\alias{sam_clsSigma,immunoMeta-method}
\alias{sam_clsEvents.immunoMeta}
\alias{sam_clsEvents,immunoMeta-method}
\alias{nobs.immunoMeta}
\alias{nobs,immunoMeta-method}
\alias{ncls.immunoMeta}
\alias{ncls,immunoMeta-method}
\alias{npar.immunoMeta}
\alias{npar,immunoMeta-method}
\alias{parameters.immunoMeta}
\alias{parameters,immunoMeta-method}
\alias{parameters<-.immunoMeta}
\alias{parameters<-,immunoMeta,character-method}
\alias{label.immunoMeta}
\alias{label,immunoMeta-method}
\alias{weights.immunoMeta}
\alias{weights,immunoMeta-method}
\alias{mu.immunoMeta}
\alias{mu,immunoMeta-method}
\alias{sigma.immunoMeta}
\alias{sigma,immunoMeta-method}
\alias{events}
\alias{events.immunoMeta}
\alias{events,immunoMeta-method}
\alias{posterior.immunoMeta}
\alias{posterior,immunoMeta-method}
\alias{prop}
\alias{prop<-}
\alias{prop.immunoMeta}
\alias{prop<-.immunoMeta}
\alias{prop,immunoMeta-method}
\alias{prop<-,immunoMeta-method}
\alias{desc}
\alias{desc<-}
\alias{desc.immunoMeta}
\alias{desc<-.immunoMeta}
\alias{desc,immunoMeta-method}
\alias{desc<-,immunoMeta-method}
\alias{descFull}
\alias{descFull.immunoMeta}
\alias{descFull,immunoMeta-method}
\alias{level}
\alias{level.immunoMeta}
\alias{level,immunoMeta-method}
\alias{level<-}
\alias{level<-.immunoMeta}
\alias{findLevel}
\alias{findLevel.immunoMeta}
\alias{findLevel,immunoMeta-method}
\alias{clusters}
\alias{classified}
\alias{unclassified}
\alias{clusters.immunoMeta}
\alias{classified.immunoMeta}
\alias{unclassified.immunoMeta}
\alias{clusters,immunoMeta-method}
\alias{classified,immunoMeta-method}
\alias{unclassified,immunoMeta-method}
\alias{addLevel<-}
\alias{addLevel<-.immunoMeta}
\alias{move<-}
\alias{move<-.immunoMeta}
\alias{parent<-}
\alias{parent<-.immunoMeta}
\alias{parent<-.immunoMeta-method}
\alias{remove<-}
\alias{remove<-.immunoMeta}
\alias{transfer<-}
\alias{transfer<-.immunoMeta}
\alias{transfer<-,immunoMeta,immunoMeta-method}
\alias{finalize}
\alias{finalize.immunoMeta}
\alias{finalize,immunoMeta-method}
\alias{subset.immunoMeta}
\alias{subset,immunoMeta-method}
\alias{transformParams.immunoMeta}
\alias{transformParams,immunoMeta-method}
\alias{clusterCoeff,immunoMeta-method}
\alias{clusterDist,immunoMeta-method}
\alias{clusterProb,immunoMeta-method}

\encoding{utf8}

\title{
Acessors and Methods for immunoMeta Objects
}

\description{
Documentation of the accessors and methods for immunoMeta-objects
}


\arguments{
\item{object, immunoMeta}{an object of class
\code{\link[=immunoMeta-object]{immunoMeta}} as return by
\code{\link{meta.process}}.}

\item{cls}{ cluster subset for retrieved slot values. }
\item{par}{ parameter subset for retrieved slot values. }

\item{pos}{Gives the position in the immunoMeta-hierarchy. \code{pos}
is an array of indices which addresses the level of interest. Each level
in the immunoMeta-hierarchy consists of a name (\code{desc}), meta-cluster
subset (array of cluster indices) and a vector of sub-levels. \code{pos} is
the sequence of indices into these sub-levels begining at root level.} 

}


\section{Accessors}{
\describe{

\item{nsam}{
    the number of \code{immunoClust}-objects (samples) which
    are co-clustered.

    \emph{Usage:}

    \code{nsam(immunoMeta)}
}

\item{sam_ncls}{
    the number of cell event clusters in the\code{immunoClust}-objects 
    (samples) which are co-clustered.

    \emph{Usage:}

    \code{sam_ncls(immunoMeta, for.samples=seq_len(nsam(meta))}
}

\item{sam_clsWeights}{
    the weigths of all cell event clusters which are collected
    for co-clustering.

    \emph{Usage:}

    \code{sam_clsWeights(immunoMeta)}
}

\item{sam_clsMu}{
    the means of all cell event clusters which are collected
    for co-clustering.

    \emph{Usage:}

    \code{sam_clsMu(immunoMeta)}
}

\item{sam_clsSigma}{
    the co-variance matrices of all cell event clusters which are collected
    for co-clustering.

    \emph{Usage:}

    \code{sam_clsSigma(immunoMeta)}
}

\item{sam_clsEvents}{
    the event numbers of all cell event clusters which are collected
    for co-clustering.

    \emph{Usage:}

    \code{sam_clsEvents(immunoMeta)}
}


\item{nobj}{
    the number of cell events clusters from sample cell-clustering which
    are co-clustered.

    \emph{Usage:}

    \code{nobj(immunoMeta)}
}

\item{ncls}{
    the number of meta-clusters.

    \emph{Usage:}

    \code{ncls(immunoMeta)}
}

\item{npar}{
    the number of parameters measured, cell-clustered and meta-clustered

    \emph{Usage:}

    \code{npar(immunoMeta)}
}

\item{parameters, parameters<-}{
    extracts or replaces the names of measured, cell-clustered and
    meta-clustered parameters

    \emph{Usage:}

    \code{parameters(immunoMeta)}

    \code{parameters(immunoMeta) <- value}
}

\item{label}{
    the meta-clustering label, that is the assignment of the cell-clusters
    to the meta-clusters. 


    \emph{Usage:}

    \code{label(immunoMeta, for.sample=NA)}

    If for.sample is specified, the label part for this sample only.
}

\item{weights}{
    the meta-clustering weights for the cluster selection (all meta-cluster by
    default)

    \emph{Usage:}

    \code{weights(immunoMets,cls=seq_len(ncls(immunoMeta)))}

}

\item{mu}{
    the meta-cluster mean values for the cluster and parameter selection 
    (all meta-cluster and all parameter by default)

    \emph{Usage:}

    \code{mu(immunoMeta, cls=seq_len(ncls(immunoMeta)), 
        par=seq_len(npar(immunoMeta)))}

}

\item{sigma}{
    the meta-cluster co-variance values for the cluster and parameter selection 
    (all meta-cluster and all parameter by default)

    \emph{Usage:}

    \code{sigma(immunoMeta, cls=seq_len(ncls(immunoMeta)),
        par=seq_len(npar(immunoMeta)))}

}

\item{posterior}{
    the posterior probabilities of cluster membership for each cell-cluster

    \emph{Usage:}

    \code{posterior(immunoMeta)}

}


\item{events}{
    the cell-event numbers for each sample for the cluster selection
    (all meta-cluster by default)

    \emph{Usage:}

    \code{events(immunoMeta, ncls=seq_len(ncls(immunoMeta)))}
}

\item{prop, prop<-}{
    get or a property value in the hierarchy level given by \code{pos} and
    named \code{name}

    \emph{Usage:}

    \code{prop(immunoMeta, name, pos=c()) }

    \code{prop(immunoMeta, name, pos, for.level=TRUE, for.sublevels=FALSE)
        <- value }

If the option \code{for.sublevels} is set, the property value will by setted 
deep for all sub-levels of the by \code{pos} specified level.

The \code{prop} interface is very basic and no checks for meaningfull 
properties and values are performed. It could be used for everything at 
any time. Nevertheless, there are some property keys which are used internally
mainly to control the plot routine for the levels.

    \code{desc}     the name of this level.

    \code{M}        the mean of all clusters in this level

    \code{S}        the co-variance meatirx of all clusters in this level

    \code{pscales}  a list of \code{npar} entries for the limits and ticksi
    information. Normaly, only set on root-level and then used for all
    sub-levels. But could set and altered at any level.

    \code{plot.subset}  an array of parameter indices used as default for the
    plot this level.

    \code{plot.color}   an index in the \code{palette} or other specified
    color used for plots of this level in its parent level.

    \code{plot.childs}  to be renamed in \code{plot.levels}.

    \code{plot.parent}  when set, additionally all cluster of the parent level
    are plotted in light gray.

}

\item{desc, desc<-}{
    Get or set the \code{desc} property in the by \code{pos} specified level.

    \emph{Usage:}

    \code{desc(immunoMeta, pos)}

    \code{desc(immunoMeta, pos) <- value}
}

\item{descFull}{
Gives the full description path for the level given by \code{pos}, i.e. the
concatinate \code{desc} values of this all parent levels.

    \emph{Usage:}

    \code{descFull(immunoMeta, pos)}

}

\item{level, level<- }{
Get or replace the level object at specified \code{pos},

    \emph{Usage:}

    \code{value <- level(immunoMeta, pos)}

    \code{level(immunoMeta, pos ) <- value}

}

\item{findLevel}{
Find the level \code{pos} value for a specific cluster \code{cls}

    \emph{Usage:}

    \code{pos <- findLevel(immunoMeta, cls)}

}

\item{clusters}{
Retrieves the cluster subset for the level at \code{pos}.

    \emph{Usage:}

    \code{cls <- clusters(immunoMeta, pos)}

}

\item{classified}{
Retrieves the cluster subset for the level at \code{pos} which are classified
in sub-levels.


    \emph{Usage:}

    \code{cls <- classified(immunoMeta, pos)}

}

\item{unclassified}{
Retrieves the cluster subset for the level at \code{pos} which are not
classified in sub-levels.


    \emph{Usage:}

    \code{cls <- unclassified(immunoMeta, pos)}
}
}
}

\section{Manipulators}{
\describe{
\item{addLevel<-}{
Adds a level at a specified hierarchy position \code{pos}. A level consists
of a name (\code{desc}) and a cluster subset \code{cls}.

    \emph{Usage:}

    \code{addLevel(immunoMeta, pos, desc="new level") <- cls}

}

\item{move<-}{
Moves a cluster subset to a specific immunoMeta level. Clusters in \code{cls}
are added to parent levels if nessesary and removed from other levels.

    \emph{Usage:}

    \code{move(immunoMeta, pos) <- cls}
}

\item{remove<-}{
removes a cluster subset from a specific immunoMeta level. 

    \emph{Usage:}

    \code{remove(immunoMeta, pos) <- cls}
}

\item{parent<-}{
sets the parent for this level, or this level as parent for all
its sub-levels

    \emph{Usage:}

    \code{parent(immunoMeta, pos) <- c()}

    \code{parent(immunoMeta, pos) <- level}
}

\item{transfer<-}{
Overtakes the annotation of an immunoMeta-object to this immunoMeta-object

    \emph{Usage:}

    \code{transfer(immunoMeta) <- annotatedMeta}

}

}
}

\section{Methods}{
\describe{
\item{finalize}{
After manipulations of a \code{immunoMeta}-object \code{finalize}
restructure all levels and returns the finalized object, where the
parent relations are solved and the mean and co-variances of all levels
are build.

    \emph{Usage:}

    \code{immunoMeta <- finalize(immunoMeta)}

}

\item{subset}{
Builds the \code{immunoMeta}-object for a cluster and/or parameter subset


    \emph{Usage:}

    \code{subsetMeta <- subset(immunoMeta, cls=seq_len(ncls(meta)),
        par=seq_len(npar(meta)))}

}

\item{transformParams}{
Scales and translates the cluster means of the \code{immunoMeta}-object
in each parameter

    \emph{Usage:}

    \code{transformedMeta <- transformParams(immunoMeta, scale=c(), offset=c())}
}

\item{clusterCoeff}{
Calculates the bhattacharrya coefficients of clusters \code{cls} for a level
\code{lvl} in the \code{immunoMeta}-object

    \emph{Usage:}

    \code{ret <- clustersCoeff(immunoMeta, cls, lvl, 
        par=seq_len(npar(immunoMeta))}

}

\item{clusterDist}{
Calculates the bhattacharrya distances of clusters \code{cls} for a level
\code{lvl} in the \code{immunoMeta}-object

    \emph{Usage:}

    \code{ret <- clustersDist(immunoMeta, cls, lvl, 
        par=seq_len(npar(immunoMeta))}

}

\item{clusterProb}{
Calculates the bhattacharrya probabilities of clusters \code{cls} for a level
\code{lvl} in the \code{immunoMeta}-object

    \emph{Usage:}

    \code{ret <- clustersProb(immunoMeta, cls, lvl, 
        par=seq_len(npar(immunoMeta))}

}

}
}

\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}


\seealso{
\code{\link[=immunoMeta-class]{immunoMeta-object}}
}

\examples{
###
data(dat.meta)
npar(dat.meta)
ncls(dat.meta)
cls <- clusters(dat.meta,c(1))
move(dat.meta,c(2)) <- cls
}

\keyword{immunoMeta}
