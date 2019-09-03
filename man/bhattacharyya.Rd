\name{bhattacharyya}
\alias{bhattacharyya.dist}
\alias{bhattacharyya.coeff}
\alias{bhattacharyya.prob}

\encoding{utf8}

\title{
Bhattacharyya Distance, Coefficient and Probability
}
\description{
Calculates the Bhattacharyya Distance, Coefficient and Probability 
}
\usage{
bhattacharyya.prob(gM,gS, cM,cS, alpha=1)

bhattacharyya.dist(gM, gS, cM, cS)

bhattacharyya.coeff(gM,gS, cM,cS, alpha=1)
}

\arguments{
    \item{gM, cM}{ \code{P}-dimensional vector of cluster means}
    \item{gS, cS}{ \code{PxP}-dimensinal matrix of clusters co-variances}
    \item{alpha}{ A value between 0 and 1 used to balance the bhattacharrya 
probabilities, coefficients calculated with either the full covariance 
matrices or using only the diagonal elements of it. }
}
\details{
Calculates the bhattacharyya probabilty, distance or coefficient of the
clusters, i.e. Gaussian distributions. Distance and Coefficient are symetric 
for both clusters, whereas the probabity is not.
}
\value{
The Bhattacharyya probability, distance or coefficient
}

\author{
Till Sörensen \email{till-antoni.soerensen@charite.de}
}
