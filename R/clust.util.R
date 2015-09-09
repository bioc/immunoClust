###
## clust.hclass
##  retrieve group membership from hcPairs structure
###
.clust.partconv <- function(x, consec = TRUE)
{
    n <- length(x)
    y <- numeric(n)
    u <- unique(x)
    if(consec) {
# number groups in order of first row appearance
        l <- length(u)
        for(i in 1:l)
        y[x == u[i]] <- i
    }
    else {
# represent each group by its lowest-numbered member
        for(i in u) {
            l <- x == i
            y[l] <- (1:n)[l][1]
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
    return(matrix(initial, ncol = 1, dimnames = list(NULL, as.character(G))))
    
    bad <- select < 0 | select >= k
    if (any(bad)) 
    stop("Some specified number of clusters are inconsistent")
    
## number of groupings to be labeled
    L <- length(select)
    cl <- matrix(NA, nrow = n, ncol = L, dimnames = list(NULL, as.character(G)))
    if (select[1]) { 
        m <- 1
    }
    else {
## highest G == n
        cl[, 1] <- initial
        m <- 2
    }
    for (l in 1:max(select)) {
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
    apply(cl[, L:1, drop = FALSE], 2, .clust.partconv, consec = TRUE)
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
#         order = integer(n+1), iia = integer(n+1), iib = integer(n+1),
#         PACKAGE="stats")
# hcls <- list(merge = cbind(obj$iia[1L:n], obj$iib[1L:n]), 
#        height = crit, crit=crit, order = obj$order, 
#        labels = 1:n)
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
    for( n in 1:(N-1) ) {
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
.immunoClust2 <- function(obj, K, P, N, expName="", parameters=c(), inc=c())
{
    L <- obj$L
    
#output sigma    
    sigma <- array(0, c(L, P, P))
    s <- matrix(obj$s, K, P * P, byrow=TRUE)
    for (k in 1:L)
    sigma[k,,] <- matrix(s[k,], P, P, byrow = TRUE)

    
# output mu    
    mu <- matrix(obj$m, K, P, byrow=TRUE)[1:L,]
    dim(mu) <- c(L,P)
    
    
# output z
    if( length(inc) > 0 ) {
        z <- matrix(NA, length(inc), K)
        z[inc,] <- matrix(obj$z, N, K, byrow=TRUE)
        z <- z[,1:L]
        
        label <- rep(NA, length(inc))
        label[inc] <- obj$label
        
        N <- length(inc)
        
    }
    else {
        z <- matrix(obj$z, N, K, byrow=TRUE)
        z <- as.matrix(z[,1:L])
## 2014.07.08: quick fix (na/nan should not occur)  
        if( sum(is.infinite(z) | is.na(z) | is.nan(z)) > 0 ) {
            warning("Fehler: Z has infinite values", 
                sum(is.infinite(z) | is.na(z) | is.nan(z) ), "\n")
            z[is.infinite(z) | is.na(z) | is.nan(z)] <- 0
        }
        
        label <- obj$label
    }
    
    dim(z) <- c(N,L)
    
    
#    if( !is.null(history) ) {
##    for( k in 1:L ) {
##      history[k] <- history[ obj$history[k] ]
##    }
#        history <- paste(sep="", history[1:L], ",", 1:L)
#    }
#    else {
#        history <- paste(1:L)
#    }
#    if( !is.null(state) ) {
##    for( k in 1:L ) {
##      state[k] <- state[ obj$history[k] ]
##    }
#        state <- state[1:L]
#    }
#    else {
#        state <- rep(0, L)
#    }
    
    new("immunoClust", expName=expName, parameters=parameters,
        K=L, P=P, N=N, w=obj$w[1:L], mu=mu, sigma=sigma,
        z=z, label=label, 
        logLike=obj$logLike, BIC=obj$logLike[1], ICL=obj$logLike[2],
        history=paste(obj$history), state=rep(0,L))
    
}


##########
## meta-cluster separation trails and utils
##########

## calculate discrimination score for X given a separation a
.separation.score <- function( X, a)
{
    names(X) <- NULL
    M <- mean(X)
    K <- length(a)
    M_g <- rep(0, K+1)
    M_g[1] <- mean(X[X <= a[1]])
    if( K > 1 ) {
        for( g in 2:K ) {
            M_g[g] <- mean(X[a[g-1] < X & X < a[g]])
        }
    }
    M_g[K+1] <- mean(X[X >= a[K]])
    
## between
    S_b <- 0
    S_w <- 0 
    
    for( i in 1:length(X) ) {
        x <- X[i]
        if( x <= a[1] ) {
            S_w <- S_w + (x-M_g[1])^2
            S_b <- S_b + (M-M_g[1])^2
        }
        
        if( x >= a[K] ) {
            S_w <- S_w + (x-M_g[K+1])^2
            S_b <- S_b + (M-M_g[K+1])^2
        }
        
        if( K > 1 ) {
            for( g in 2:K ) {
                if( a[g-1] < x && x < a[g] ) {
                    S_w <- S_w + (x-M_g[g])^2
                    S_b <- S_b + (M-M_g[g])^2
                }
            }
        }
    }
    
    r <- S_b / S_w
    r
}


## calculate discrimination score for X with square sigma given a separation a
.classification.score <- function( X, sigma, a)
{
    names(X) <- NULL
    M <- mean(X)
    S <- sum(sigma)
    K <- length(a)
    
    r <- 1
    if( K > 0 ) {
        M_g <- rep(0, K+1)
        M_g[1] <- mean(X[X <= a[1]])
        if( K > 1 ) {
            for( g in 2:K ) {
                M_g[g] <- mean(X[a[g-1] < X & X < a[g]])
            }
        }
        M_g[K+1] <- mean(X[X >= a[K]])
        
        S_b <- 0
        S_w <- 0 
        
        for( i in 1:length(X) ) {
            x <- X[i]
            if( x <= a[1] ) {
                S_w <- S_w + (x-M_g[1])^2
                S_b <- S_b + (M-M_g[1])^2
            }
            
            if( x >= a[K] ) {
                S_w <- S_w + (x-M_g[K+1])^2
                S_b <- S_b + (M-M_g[K+1])^2
            }
            
            if( K > 1 ) {
                for( g in 2:K ) {
                    if( a[g-1] < x && x < a[g] ) {
                        S_w <- S_w + (x-M_g[g])^2
                        S_b <- S_b + (M-M_g[g])^2
                    }
                }
            }
        }
        
        r <- 1 + S_b / (S + S_w)
    }
    
    r
}


.separateX <- function(X, t)
{
    l <- min(X)
    r <- max(X)
    for( x in X ) {
        if( x <= t && x>l )
        l <- x
        if( x >= t && x<r )
        r <- x
    }
    
    (l+r)/2
}

.minDelta <- function(X)
{
    if( length(X) <= 1 ) {
        delta <- 0
    }
    else {
        delta <- abs(X[1] - X[2])
        
        for( x in X ) {
            for( y in X ) {
                if( x != y ) {
                    if( abs(x-y) < delta ) {
                        delta <- abs(x-y)
                    }
                }
            }
        }
    }
    delta/2
    
}

.separation.max <- function( 
X, step=10000, a=c(min(X),max(X)), 
a.min=rep(min(X),length(a)), a.max=rep(max(X),length(a)) 
) {
    
    t <- a
    n <- length(a)
    
## cat("step ", step, " N=", n, "X=", length(X), "\n" )
    
## upper bound  
    if( n > 1 ) {
        t[n] <- a.max[n]
        
        N <- max(1,round((a.max[n-1] - a.min[n-1])/ step))
        S <- rep(.separation.score(X,t), N)
        
        t[n] <- max(t[n-1],t[n] - step)
        i <- 2
        while( t[n] > max(t[n-1],a.min[n]) ) {
            S[i] <- .separation.score(X, t)
            i <- i+1
            t[n] <- t[n] - step
        }
        i <- which.max(S)
        t[n] <- .separateX(X, a.max[n] - (i-1)*step)
        
        
## lower bound  
        t[1] <- a.min[1]
        N <- max(1,round((a.max[1] - a.min[1])/ step))
        S <- rep(.separation.score(X,t), N)
        
        
        t[1] <- t[1] + step
        i <- 2
        while( t[1] < min(a.max[1],t[2]) ) {
            S[i] <- .separation.score(X, t)
            i <- i+1
            t[1] <- t[1] + step
        }
        i <- which.max(S)
        t[1] <- .separateX(X, a.min[1] + (i-1)*step)
    }
    else {
## n == 1
        t <- a.min
        
        N <- max(1,round((a.max - a.min)/ step))
        S <- rep(.separation.score(X,t), N)
        
        t <- t + step
        i <- 2
        while( t <= a.max ) {
            S[i] <- .separation.score(X, t)
            i <- i+1
            t <- t + step
        }
        
        i <- which.max(S)
        t <- .separateX(X, a.min + (i-1)*step)
        
    }
## between
    if( n > 2 ) {
        for( j in 2:(n-1) ) {
            N <- max(1,round((a.max[j]-a.min[j])/step))
            if( N > 1 ) {
                a.min[j] <- max(a.min[j],t[j-1])
                t[j] <- a.min[j]
                
                S <- rep(.separation.score(X,t),N)
                
                t[j] <- t[j] + step
                i <- 2
                while( t[j] < min(t[j+1],a.max[j]) ) {
                    S[i] <- .separation.score(X, t)
                    i <- i+1
                    t[j] <- t[j] + step
                }
                i <- which.max(S)
                t[j] <- .separateX(X, a.min[j] + (i-1)*step)
                
            }
        }
    }
    
    t
}

.separation.maximization <- 
function(X, step.size=10000, step.min=1, more=FALSE)
{
    
    delta <- .minDelta(X)
    
    a.min <- min(X)+delta
    a.max <- max(X)-delta
    
    a <- .separation.max(X, a=c(a.min, a.max), step=step.size)
    b <- .separation.max(X, a=a.min, a.min=a.min, a.max=a.max, step=step.size)
    
    if( more ) {
        a <- .separation.max(X, a=c(min(X), a), a.min=c(min(X), a-step.size), 
                            a.max=c(max(X), a+step.size), step=step.size )
    }
    
    while(step.size > max(step.min,delta) ) {
##    cat("step", step.size, "\n")
        a <- .separation.max(X, a=a, a.min=a-step.size, 
                            a.max=a+step.size, step=step.size/10)
        b <- .separation.max(X, a=b, a.min=max(a.min,b-step.size), 
                            a.max=min(a.max,b+step.size), step=step.size/10)
        step.size <- step.size/10
    }
    
    G <- length(a)
    if( a[1] == min(X) ) {
        a[1] <- min(X) + delta
    }
    if( a[G] == max(X) ) {
        a[G] <- max(X) - delta
    }
    
    if( G == 2 ) {
        aScore <- .separation.score(X,a)
        bScore <- .separation.score(X,b)
        
        if( aScore==0 ) {
            cat("\t", X, "\n")
        }
        
        if( 2*bScore >= aScore ) {
## last gate has no impact: -> no dim
            a[1] <- a[2] <- b 
        }
        
    }
    
    a
    
}

.classification.max <- function( 
X, S, step=10000, a=c(min(X),max(X)), 
a.min=rep(min(X),length(a)), a.max=rep(max(X),length(a)) 
) {
    
    t <- a
    n <- length(a)
    
    
## upper bound  
    if( n > 1 ) {
        t[n] <- a.max[n]
        
        N <- max(1,round((a.max[n-1] - a.min[n-1])/ step))
        score <- rep(.classification.score(X,S,t), N)
        
        t[n] <- max(t[n-1],t[n] - step)
        i <- 2
        while( t[n] > max(t[n-1],a.min[n]) ) {
            score[i] <- .classification.score(X, S, t)
            i <- i+1
            t[n] <- t[n] - step
        }

        i <- which.max(score)
        t[n] <- .separateX(X, a.max[n] - (i-1)*step)
        
## lower bound  
        t[1] <- a.min[1]
        N <- max(1,round((a.max[1] - a.min[1])/ step))
        score <- rep(.classification.score(X,S,t), N)
        
        t[1] <- t[1] + step
        i <- 2
        while( t[1] < min(a.max[1],t[2]) ) {
            score[i] <- .classification.score(X, S, t)
            i <- i+1
            t[1] <- t[1] + step
        }

        i <- which.max(score)
        t[1] <- .separateX(X, a.min[1] + (i-1)*step)
    }
    else {
## n == 1
        t <- a.min
        
        N <- max(1,round((a.max - a.min)/ step))
        score <- rep(.classification.score(X,S,t), N)
        
        t <- t + step
        i <- 2
        while( t <= a.max ) {
            score[i] <- .classification.score(X, S, t)
            i <- i+1
            t <- t + step
        }
        
        i <- which.max(score)
        t <- .separateX(X, a.min + (i-1)*step)
        
    }
## between
    if( n > 2 ) {
        for( j in 2:(n-1) ) {
            N <- max(1,round((a.max[j]-a.min[j])/step))
            if( N > 1 ) {
                a.min[j] <- max(a.min[j],t[j-1])
                t[j] <- a.min[j]
                
                score <- rep(.classification.score(X,S,t),N)
                
                t[j] <- t[j] + step
                i <- 2
                while( t[j] < min(t[j+1],a.max[j]) ) {
                    score[i] <- .classification.score(X, S, t)
                    i <- i+1
                    t[j] <- t[j] + step
                }
                i <- which.max(score)
                t[j] <- .separateX(X, a.min[j] + (i-1)*step)
                
            }
        }
    }
    
    t
}

.classification.maximization <- function(
X, S, step.size=1, step.min=1e-5, aFac=7/6, bFac=7/6.25
) {
    
    delta <- .minDelta(X)
    
    a.min <- min(X)+delta
    a.max <- max(X)-delta
    
    a <- .classification.max(X, S, a=c(a.min, a.max), step=step.size)
    b <- .classification.max(X, S, a=(a.min+a.max)/2, 
                            a.min=a.min, a.max=a.max, step=step.size)
    
    
    while(step.size > step.min ) {
        a <- .classification.max(X, S, a=a, a.min=a-step.size, 
                                a.max=a+step.size, step=step.size/10)
        b <- .classification.max(X, S, a=b, a.min=max(a.min,b-step.size), 
                                a.max=min(a.max,b+step.size), 
                                step=step.size/10)
        step.size <- step.size/10
    }
    
    G <- length(a)
    if( a[1] == min(X) ) {
        a[1] <- min(X) + delta
    }
    if( a[G] == max(X) ) {
        a[G] <- max(X) - delta
    }
    
    aScore <- .classification.score(X,S,a)
    bScore <- .classification.score(X,S,b)
    
    if( aScore <= aFac*bScore ) {
## last gate has no impact: -> no dim
        a[1] <- a[2] <- b 
    }
    if( bScore <= bFac ) {
## classification has no inpact: ->no gates 
        a[1] <- a[2] <- NA
    }
    
    a
    
}

.classification.matrices <- function(mu1, sigma1, mu2, sigma2)
{
    P <- ncol(mu1)
    
    N_1 <- nrow(mu1)
    N_2 <- nrow(mu2)
    N_g <- N_1+N_2
    M_1 <- rep(0.0, P)
    M_2 <- rep(0.0, P)
    M_g <- rep(0.0, P)

    S_b <- matrix(0.0, nrow=P, ncol=P)
    S_w <- matrix(0.0, nrow=P, ncol=P)
    
    for( i in 1:N_1 ) {
        M_1 <- M_1 + mu1[i,]
    }
    for( i in 1:N_2 ) {
        M_2 <- M_2 + mu2[i,]
    }
    
    M_g <- (M_1 + M_2) / (N_1 + N_2)
    M_1 <- M_1 / N_1
    M_2 <- M_2 / N_2
    
    for( i in 1:N_1 ) {
        X <- mu1[i,] - M_1
        S_w <- S_w + sigma1[i,,] + X %*% t(X) 
    }
    
    for( i in 1:N_2 ) {
        X <- mu2[i,] - M_2
        S_w <- S_w + sigma2[i,,] + X %*% t(X) 
    }
    S_b <- N_1 * ((M_1 - M_g) %*% t(M_1 - M_g)) + 
            N_2 * ((M_2 - M_g) %*% t(M_2 - M_g)) 
    
    e <- eigen(S_w + S_b)
    
    A <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
    e <- eigen(A %*% solve(S_w) %*% A)
    
    v <- solve(A) %*% e$vectors[,1]
    
    list("S_b"=S_b, "S_w"=S_w, "M_1"=M_1, "M_2"=M_2, "N_1"=N_1, "N_2"=N_2, 
        "A"=A, "center"=M_g, "vector"=v, "value"=e$values[1], "eigen"=e)
    
}

##
# classification separation
##
.classification.separation <- function(M, S, label )
{
    P <- ncol(M)
    m1 <- rep(0, P)
    m2 <- rep(0, P)
    s1 <- matrix(0, ncol=P, nrow=P)
    s2 <- matrix(0, ncol=P, nrow=P)
    
    for( i in which(label==1) ) {
        m1 <- m1 + M[i,]
    }
    m1 <- m1 / sum(label==1)
    for( i in which(label>1) ) {
        m2 <- m2 + M[i,]
    }
    m2 <- m2 / sum(label>1)
    for( i in which(label==1) ) {
        s1 <- s1 + S[i,,] + (M[i,]-m1) %*% t(M[i,]-m1)
    }
    s1 <- s1 /sum(label==1)
    for( i in which(label>1) ) {
        s2 <- s2 + S[i,,] + (M[i,]-m2) %*% t(M[i,]-m2)
    }
    s2 <- s2 / sum(label>1)
    d1 <- 0
    try(d1 <- mahalanobis(m1,m2, (s1 + s2)/2))
    
    m1 <- rep(0, P)
    m2 <- rep(0, P)
    s1 <- matrix(0, ncol=P, nrow=P)
    s2 <- matrix(0, ncol=P, nrow=P)
    
    for( i in which(label < 3) ) {
        m1 <- m1 + M[i,]
    }
    m1 <- m1 / sum(label < 3)
    for( i in which(label==3) ) {
        m2 <- m2 + M[i,]
    }
    m2 <- m2 / sum(label==3)
    for( i in which(label<3) ) {
        s1 <- s1 + S[i,,] + (M[i,]-m1) %*% t(M[i,]-m1)
    }
    s1 <- s1 /sum(label<3)
    for( i in which(label==3) ) {
        s2 <- s2 + S[i,,] + (M[i,]-m2) %*% t(M[i,]-m2)
    }
    s2 <- s2 / sum(label==3)
    d2 <- 0
    try(d2 <- mahalanobis(m1,m2, (s1 + s2)/2))
    
    (d1+d2)/2
}
## classification separation 


### LDA utils

