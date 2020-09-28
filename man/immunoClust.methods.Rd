\name{methods.immunoClust}
\docType{class}
\alias{nobs.immunoClust}
\alias{nobs,immunoClust-method}
\alias{ncls.immunoClust}
\alias{ncls,immunoClust-method}
\alias{npar.immunoClust}
\alias{npar,immunoClust-method}
\alias{parameters.immunoClust}
\alias{parameters,immunoClust-method}
\alias{parameters<-.immunoClust}
\alias{parameters<-,immunoClust,character-method}
\alias{label.immunoClust}
\alias{label,immunoClust-method}
\alias{weights.immunoClust}
\alias{weights,immunoClust-method}
\alias{mu.immunoClust}
\alias{mu,immunoClust-method}
\alias{sigma.immunoClust}
\alias{sigma,immunoClust-method}
\alias{events.immunoClust}
\alias{events,immunoClust-method}
\alias{posterior.immunoClust}
\alias{posterior,immunoClust-method}

\alias{subset.immunoClust}
\alias{subset,immunoClust-method}
\alias{transformParams.immunoClust}
\alias{transformParams,immunoClust-method}

\encoding{utf8}

\title{
Acessors and Methods for immunoClust Objects
}

\description{  
Documentation of the accessors and methods for immunoClust-objects
}


\arguments{
\item{object, immunoClust}{an object of class 
\code{\link[=immunoClust.object]{immunoClust}} as return by 
\code{\link{cell.process}}.}

\item{cls}{ cluster subset for retrieved slot values. }
\item{par}{ parameter subset for retrieved slot values. }

}


\section{Accessors}{
\describe{

\item{nobs}{
    the number of cell events clustered

    \emph{Usage:}

    \code{nobs(immunoClust)}
}

\item{ncls}{
    the number of clusters.

    \emph{Usage:}

    \code{ncls(immunoClust)}
}

\item{npar}{
    the number of parameters measured, cell-clustered 

    \emph{Usage:}

    \code{npar(immunoClust)}
}

\item{parameters, parameters<-}{
    extracts or replaces the names of measured, cell-clustered parameters

    \emph{Usage:}

    \code{parameters(immunoClust)}

    \code{parameters(immunoClust) <- value}
}

\item{label}{
    the clustering label, that is the assignment of the cell-events to the 
    clusters. 


    \emph{Usage:}

    \code{label(immunoClust)}


}

\item{weights}{
    the clustering weights for the cluster selection (all cluster by default)

    \emph{Usage:}

    \code{weights(immunoClust,cls=seq_len(ncls(immunoClust)))}

}

\item{mu}{
    the cluster mean values for the cluster and parameter selection 
    (all cluster and all parameter by default)

    \emph{Usage:}

    \code{mu(immunoClust, cls=seq_len(ncls(immunoClust)), 
            par=seq_len(npar(immunoClust)))}

}

\item{sigma}{
    the cluster co-variance values for the cluster and parameter selection 
    (all cluster and all parameter by default)

    \emph{Usage:}

    \code{sigma(immunoClust, cls=seq_len(ncls(immunoClust)), 
        par=seq_len(npar(immunoClust)))}

}

\item{posterior}{
    the posterior probabilities of cluster membership for each event

    \emph{Usage:}

    \code{posterior(immunoClust)}

}

\item{events}{
    the cell-event numbers for the cluster selection (all cluster by default)

    \emph{Usage:}

    \code{events(immunoClust, ncls=seq_len(ncls(immunoClust)))}
}

}
}

\section{Methods}{
\describe{
    \item{subset}{
    Builds the \code{immunoClust}-object for a parameter subset


    \emph{Usage:}

    \code{res <- subset(immunoClust, par)}

}

    \item{transformParams}{
    Scales and translates the cluster means of the \code{immunoClust}-object
    in each parameter

    \emph{Usage:}

    \code{res <- transformParams(immunoClust, scale=c(), offset=c())}
}

}
}

\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}


\seealso{
\code{\link[=immunoClust-class]{immunoClust-object}}
}

\examples{
###
data(dat.exp)
## cell.clustering result for dat.fcs
res <- dat.exp[[1]]
nobs(res)
ncls(res)

}

\keyword{immunoClust}
