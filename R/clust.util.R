
##
bhattacharyya.prob <- function(gM,gS, cM,cS, alpha=1)
{
    ret <- 0
    if( is.null(gS) || is.null(cS) )
        return (0)
        
    gM <- as.vector(gM)
    cM <- as.vector(cM)
    gS <- as.matrix(gS)
    cS <- as.matrix(cS)
        
    if( alpha <= 0 ) {
        ret <- bhattacharyya.prob(gM,diag(diag(gS)), cM, diag(diag(cS)))
        return(ret)
    }
    else if( alpha < 1 ) {
        a <- bhattacharyya.prob(gM,gS, cM, cS)
        b <- bhattacharyya.prob(gM,diag(diag(gS)), cM, diag(diag(cS)))
        ret <- alpha*a + (1-alpha)*b
        return(ret)
    }
    
    
        
    suppressWarnings(det_g <- log(det(gS)))
    suppressWarnings(det_c <- log(det(cS)))
    # use chol robustness?
    #det_g <- 2*log(det(chol(gS)))
    #det_c <- 2*log(det(chol(cS)))
    
    S <- NULL
    suppressWarnings(S <- chol((gS+cS)/2))
    
    if( is.null(S) ) {
        return(0)
    }
    
    S <- chol2inv(S)
    
    logD <- log(det(S))
    dist <- mahalanobis(gM, cM, S, inverted=TRUE)
    logD <- logD + 0.5*det_g + 0.5*det_c
    logD <- logD - 0.5*0.5*dist
    
    ## normalization factor
    logD <- logD - 0.25*det_g
    
    ret <- exp(0.5*logD)
    if( is.na(ret) || ret <= 0 ) {
        ret <- 0
    }
    
    ret
    
}

## bhattacharyya.dist = -log(bhattacharyya.coeff)
#bhattacharyya.dist <- function (gM, gS, cM, cS)
#{
#   if( is.null(gS) || is.null(cS) ) {
#
#        return(0)
#    }
#
#    S <- (gS + cS)/2
#    d1 <- mahalanobis(gM, cM, S)/8
#
#    d2 <- log(det(as.matrix(S))/sqrt(det(as.matrix(gS)) *
#    det(as.matrix(cS))))/2
#    ret <- d1 + d2
#    ret
#}

bhattacharyya.dist <- function (gM, gS, cM, cS)
{
    if( is.null(gS) || is.null(cS) ) {
        return(0)
    }
    
    gM <- as.vector(gM)
    cM <- as.vector(cM)
    gS <- as.matrix(gS)
    cS <- as.matrix(cS)
    
    ## use chol: robustness? anyway
    suppressWarnings(det_g <- log(det(gS)))
    suppressWarnings(det_c <- log(det(cS)))
    #det_g <- 2*log(det(chol(gS)))
    #det_c <- 2*log(det(chol(cS)))
    
    
    S <- NULL
    try( S <- solve((gS+cS) / 2), silent=TRUE)
    
    if( is.null(S) ) {
        return (0)
    }
    
    d1 <- mahalanobis(gM, cM, S, inverted=TRUE)/8
    d2 <- (log(det(S)) + 0.5*det_g + 0.5*det_c)/2
    
    ret <- d1 - d2
    ret
}

## bhattacharyya.coeff = exp(-bhattacharyya.dist)
bhattacharyya.coeff <- function(gM,gS, cM,cS, alpha=1)
{
    ret <- 0
    if( is.null(gS) || is.null(cS) ) {
        return(0)
    }
    
    gM <- as.vector(gM)
    cM <- as.vector(cM)
    gS <- as.matrix(gS)
    cS <- as.matrix(cS)
    
    if( alpha <= 0 ) {
        ret <- bhattacharyya.coeff(gM,diag(diag(gS)), cM, diag(diag(cS)))
        return(ret)
    }
    else if( alpha < 1 ) {
        a <- bhattacharyya.coeff(gM,gS, cM, cS)
        b <- bhattacharyya.coeff(gM,diag(diag(gS)), cM, diag(diag(cS)))
        ret <- alpha*a + (1-alpha)*b
        return(ret)
    }
    
    
    suppressWarnings(det_g <- log(det(gS)))
    suppressWarnings(det_c <- log(det(cS)))
    # use chol robuster? anyway something strange if det < 0
    #det_g <- 2*log(det(chol(gS)))
    #det_c <- 2*log(det(chol(cS)))
    
    S <- NULL
    try( S <- solve((gS+cS) / 2), silent=TRUE)
    
    if( is.null(S) ) {
        return (0)
    }
    
    logD <- log(det(S))
    dist <- mahalanobis(gM, cM, S, inverted=TRUE)
    logD <- logD + 0.5*det_g + 0.5*det_c
    logD <- logD - 0.5*0.5*dist
    
    ret <- exp(0.5*logD)
    if( is.na(ret) || ret <= 0 ) {
        ret <- 0
    }
    
    if( ret > 1 ) {
        ret <- 1
    }
    
    ret
    
}
####

###
## clust.hclass
##  retrieve group membership from hcPairs structure
###
.clust.partconv <- function(x, consec=TRUE)
{
    n <- length(x)
    y <- numeric(n)
    u <- unique(x)
    if(consec) {
# number groups in order of first row appearance
        l <- length(u)
        for(i in seq_len(l))
        y[x == u[i]] <- i
    }
    else {
# represent each group by its lowest-numbered member
        for(i in u) {
            l <- x == i
            #y[l] <- (1:n)[l][1]
            y[l] <- seq_len(n)[l][1]
        }
    }
    y
}
.clust.hclass <- function (hcPairs, G) 
{
    initial <- attributes(hcPairs)$init
    n <- length(initial)
    k <- length(unique(initial))
    G <- if (missing(G)) k:2 else rev(sort(unique(G)))
    select <- k - G
    
    if (length(select) == 1 && !select) 
    return(matrix(initial, ncol=1, dimnames=list(NULL, as.character(G))))
    
    bad <- select < 0 | select >= k
    if (any(bad)) 
    stop("Some specified number of clusters are inconsistent")
    
## number of groupings to be labeled
    L <- length(select)
    cl <- matrix(NA, nrow=n, ncol=L, dimnames=list(NULL, as.character(G)))
    if (select[1]) { 
        m <- 1
    }
    else {
## highest G == n
        cl[, 1] <- initial
        m <- 2
    }
    for ( l in seq_len(max(select)) ) {
# merge at l    
        ij <- hcPairs[, l]
        i <- ij[1]
        j <- ij[2]
# i < j: all j became i
        initial[initial == j] <- i
        if (select[m] == l) {
            cl[, m] <- initial
            m <- m + 1
        }
    }
    apply(cl[, L:1, drop=FALSE], 2, .clust.partconv, consec=TRUE)
}
## clust.hclass

###
## convert hcPairs structure to hclust class
#.clust.hclust2 <- function(hcPairs)
#{
# li <- hcPairs[1,]
# lj <- hcPairs[2,]
# crit <- attributes(hcPairs)$change
# n <- length(li)
# obj <- .Fortran("hcass2", n=as.integer(n+1), 
#         ia=as.integer(li), ib=as.integer(lj),
#         order=integer(n+1), iia=integer(n+1), iib=integer(n+1),
#         PACKAGE="stats")
# hcls <- list(merge=cbind(obj$iia[1L:n], obj$iib[1L:n]), 
#        height=crit, crit=crit, order=obj$order, 
#        labels=1:n)
# class(hcls) <- "hclust"
# hcls
#}
###

###
## convert hclust merging nodes to hcPairs merging nodes 
.clust.hpairs2 <- function(li, lj)
{
    hcp <- cbind(li,lj)
    N <- nrow(hcp)
    if( N>1)
    for( n in seq_len(N-1) ) {
        i <- hcp[n,1]
        j <- hcp[n,2]
        if( i > j ) {
            hcp[n,1] <- j
            hcp[n,2] <- i
            k <- j
        }
        else {
            k <- i
        }
        
        for( l in (n+1):N ) {
            if( hcp[l,1] == -n ) {
                hcp[l,1] <- k 
            }
            if( hcp[l,2] == -n ) {
                hcp[l,2] <- k
            }
        }
        
    }
    hcp
}
###

###
##


.clust.mergedClusters <- function(x, cls)
{
    P <- ncol(x@mu)
    M <- rep(0, P)
    S <- matrix(0, nrow=P, ncol=P)
    
    if( length(cls) > 1 ) {
        M <- colSums(x@w[cls] * x@mu[cls,]) / sum(x@w[cls])
    }
    else {
        M <- x@mu[cls,]
    }
    
    for( i in cls ) {
        for( p in seq_len(P) ) {
            for( q in seq_len(P) ) {
                S[p,q] <- S[p,q] + (x@w[i]) * ( x@sigma[i, p, q] + 
                                    (x@mu[i,p]-M[p])*(x@mu[i,q]-M[q]) ) 
            }
        }
    }    
    S <- S/sum(x@w[cls])
    
    list("mu"=M,"sigma"=S)
}

## min KL merged cluster
.clust.mergedClusters2 <- function(w, mu_, sigma_, cls)
{
    if( length(cls) == 0 ) {
        return( list("mu"=NULL, "sigma"=NULL) )
    }
    
    P <- ncol(mu_)
    M <- rep(0, P)
    S <- matrix(0, nrow=P, ncol=P)
    
    if( length(cls) > 1 ) {
        for( p in seq_len(P) )
        M[p] <- sum(w[cls] * mu_[cls,p]) / sum(w[cls])
    }
    else {
        M <- mu_[cls,]
    }
    
    for( i in cls ) {
        for( p in seq_len(P) ) {
            for( q in seq_len(P) ) {
                S[p,q] <- S[p,q] + (w[i]) * ( sigma_[i, p, q] +
                                    (mu_[i,p]-M[p])*(mu_[i,q]-M[q]) )
            }
        }
    }
    S <- S/sum(w[cls])
    
    list("mu"=M,"sigma"=S, "weight"=sum(w[cls]))
}
## min KL merged cluster

###
##
.immunoClust2 <- function(obj, K, P, N, state=NULL,
expName="", parameters=c(), inc=c())
{
    L <- obj$L
    
#output sigma    
    sigma <- array(0, c(L, P, P))
    s <- matrix(obj$s, K, P * P, byrow=TRUE)
    for (k in seq_len(L))
    sigma[k,,] <- matrix(s[k,], P, P, byrow=TRUE)

    
# output mu    
    mu <- matrix(obj$m, K, P, byrow=TRUE)[seq_len(L),]
    dim(mu) <- c(L,P)
    
    
# output z
    if( length(inc) > 0 ) {
        z <- matrix(NA, length(inc), K)
        z[inc,] <- matrix(obj$z, N, K, byrow=TRUE)
        z <- z[,seq_len(L)]
        
        label <- rep(NA, length(inc))
        label[inc] <- obj$label
        
        N <- length(inc)
        
    }
    else {
        z <- matrix(obj$z, N, K, byrow=TRUE)
        z <- as.matrix(z[,seq_len(L)])
## 2014.07.08: quick fix (na/nan should not occur)  
        if( sum(is.infinite(z) | is.na(z) | is.nan(z)) > 0 ) {
            warning("Fehler: Z has infinite values", 
                sum(is.infinite(z) | is.na(z) | is.nan(z) ), "\n")
            z[is.infinite(z) | is.na(z) | is.nan(z)] <- 0
        }
        
        label <- obj$label
    }
    
    dim(z) <- c(N,L)
    
    
    if( !is.null(state) ) {
        for( k in seq_len(L) ) {
            state[k] <- state[ obj$history[k] ]
        }
        state <- state[seq_len(L)]
    }
    else {
        state <- rep(0, L)
    }
    
    ret <-new("immunoClust", expName=expName, parameters=parameters,
        K=L, P=P, N=N, w=obj$w[seq_len(L)], mu=mu, sigma=sigma,
        z=z, label=label, 
        logLike=obj$logLike, BIC=obj$logLike[1], ICL=obj$logLike[2],
        history=paste(obj$history), state=state)
        
    attr(ret,"iterations") <- obj$iterations
    attr(ret,"tolerance") <- obj$tolerance
    
    ret
    
}
## .immunoClust2

## Default_Scales
Default_Scales <- function(trans.a, limits)
{
    if( is.null(trans.a) || is.null(limits) )
        return (NULL)
    
    pscal <- list(length(trans.a))
    for( p in seq_len(length(trans.a))) {
        if( trans.a[p] == 0 ) {
            decades <- round(log(limits[2,p],10))
            scale <- 10^(decades-2)
            minor <- 50
            ticks <- round(limits[2,p]/scale/minor)
            if( ticks < 4 ) {
                minor <- 25
                ticks <- round(limits[2,p]/scale/minor)
            }
            at <- c(0,seq_len(ticks))*minor
            labels <- sprintf("%s", at)
            i<-2
            while(i<=length(labels)) {
                labels[i] <- ""
                i <- i+2
            }
            pscal[[p]] <- list(
            at=at*scale,
            labels=labels,
            limits=limits[,p],
            unit=sprintf("[/%d]", scale)
            )
        }
        else {
            a <- trans.a[p]
            decades <- round(log(limits[2,p],10))
            linear <- abs(round(log(trans.a[p],10)))
            if( decades < linear )
                return (NULL)
            at <- c(asinh(-10^linear * a), 0, asinh(10^linear * a) )
            label <- c("", "0", "")
            for( d in (seq_len(decades-linear)+(linear-1)) ) {
                at <- c(at, asinh(5*10^d * a), asinh(10^(d+1) * a))
                label <- c(label, "", parse(text=sprintf("10^%d", d+1)) )
            }
            
            pscal[[p]] <- list(
                at= at,
                labels= label,
                limits=c(asinh(-10^linear * a),asinh(limits[2,p]*a))
                )
            
        }
    }
    pscal
    
}
## Default_Scales



