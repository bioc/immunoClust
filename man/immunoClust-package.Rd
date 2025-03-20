\name{immunoClust-package}
\alias{immunoClust-package}
\alias{immunoClust.package}
\alias{immunoClust}
\docType{package}

\encoding{utf8}

\title{
immunoClust - Automated Pipeline for Population Detection in Flow Cytometry
}
\description{

Model based clustering and meta-custering routines for Flow Cytometry (FC) data.

The immunoClust-pipeline consits of two  major procedures:

\tabular{ll}{
\code{\link{cell.process}} \tab Clustering of cell-events\cr 
\code{\link{meta.process}} \tab Meta-clustering of cell-clusters\cr 
}

Cell-events clustering is performed for each FC data sample separately. 
After this all cell-clustering results are collected in a vector and 
meta-clustering is performed to obtain the across samples popluations.
}

\author{
Till Sörensen <till-antoni.soerensen@charited.de>  
}

\references{
Sörensen, T., Baumgart, S., Durek, P., Grützkau, A. and Häupl, T. 
immunoClust - an automated analysis pipeline for the identification of 
immunophenotypic signatures in high-dimensional cytometric datasets. 
\emph{Cytometry A} (accepted). 
}

\keyword{ package }

