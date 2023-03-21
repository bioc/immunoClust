\name{dat.fcs}
\alias{dat.fcs}
\docType{data}
\title{
immunoClust Cell-clustering Sample
}
\description{
\code{flowFrame} data sample with 10.000 events in 7 parameters.
}
\usage{data(dat.fcs)}

\value{
A \code{flowCore} \code{flowFrame} with 10.000 observations on the following 
7 parameters.
\describe{
\item{\code{FCS-A}}{Forward scatter}
\item{\code{SSC-A}}{Sideward scatter}
\item{\code{FITC-A}}{CD14}
\item{\code{PE-A}}{CD19}
\item{\code{APC-A}}{CD15}
\item{\code{APC-Cy7-A}}{CD4}
\item{\code{Pacific Blue-A}}{CD3}
}
}
\details{
This FCS sample is a reduced (10.000 events) dataset in flowFrame format of the 
first sample in the dataset of immunoClust, MACS-depleted populations datasets 
2010. URL http://flowrepository.org/id/FR-FCM-ZZWB.
}
\source{
http://flowrepository.org/id/FR-FCM-ZZWB
}

\examples{
data(dat.fcs)
show(dat.fcs)
\dontrun{
## process cell clustering
dat.res <- cell.process(dat.fcs)

## apply asinh-transformation
dat.fcs.transformed <- trans.ApplyToData(dat.res, dat.fcs)

## plot results
splom(dat.res, dat.fcs.transformed,N=1000)
}

}

\keyword{datasets}
