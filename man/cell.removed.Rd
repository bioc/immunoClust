\name{cell.removed}
\alias{removed.above}
\alias{removed.below}

\encoding{utf8}

\title{
Brief Information of removed Cell-events by immunoClust Cell-event Clustering
}

\description{
Gives information about the amount of overexposed cell-event observetion in 
a FCS-file.
}

\usage{
removed.above(fcs, parameters=NULL, N=NULL, max.count=10, max=NULL)
removed.below(fcs, parameters=NULL, N=NULL, min.count=10, min=NULL)
}


\arguments{
\item{fcs}{An object of class flowFrame. Rows correspond to observations and 
columns correspond to measured parameters.}
\item{parameters}{A character vector specifying the parameters (columns) to be 
included in clustering. When it is left unspecified, all the parameters will be
used.}
\item{N}{Maximum number of observations used for clustering. When unspecified 
or higher than the number of observations (i.e. rows) in dat, all observations 
are used for clustering, otherwise only the first \code{N} observations.}
\item{max.count}{An integer specifying the threshold count for filtering data 
points from above. The default is 10, meaning that if 10 or more data points 
are larger than or equal to \code{max}, they will be excluded from the 
analysis. If \code{max} is \code{NULL}, then the maximum value of each 
parameter will be used. To suppress filtering, it is set to -1.}
\item{max}{The upper limit set for data filtering. Note that it is a vector of 
length equal to the number of parameters (columns), implying that a different 
value can be set for each parameter.}
\item{min.count}{analoguous to \code{max.count}.}
\item{min}{analoguous to \code{min}.}
}

\value{
A table with two rows containing the number of events above \code{max} in 
each parameter and above in only this parameter. The two last columns 
give the sum and percentage of all events above \code{max} in any parameter.
}

\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}

\examples{
data(dat.fcs)
removed.above(dat.fcs)
}
