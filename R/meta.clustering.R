####
# meta.process
####
meta.process <- function(
exp, dat.subset=c(), meta.iter=10, tol=1e-5, meta.bias=0.2, 
meta.alpha=.5, norm.method=0, norm.blur=2, norm.minG=10,
scatter.subset=c(), scatter.bias=0.25,scatter.prior=6
) {
    dat <- meta.exprs(exp, sub=dat.subset)
    
    res <- meta.Clustering(dat$P, dat$N, dat$K, 
                        dat$clsEvents, dat$M, dat$S, 
                        bias=meta.bias, I.iter=meta.iter, B=50, tol=tol, 
                        EM.method=20, alpha=meta.alpha, 
                        norm.method=norm.method, norm.blur=norm.blur, 
                        norm.minG=norm.minG)

    dat.norm <- dat
    if( norm.method > 0 ) {
        dat.norm <- meta.Normalize(dat$P, dat$N, dat$K, dat$clsEvents, 
                                dat$M, dat$S, res@K, res@z, 
                                method=norm.method)
    }
    
    attr(res, "desc") <- dat$desc
    
    if( length(scatter.subset) > 0 ) {
        scatter <- .meta.scatter.gating(exp, sub=scatter.subset, 
                                    EM.method=10, 
                                    EM.bias=scatter.bias, 
                                    scatter.prior=scatter.prior)
        
        res.scatter <- scatter$res
        dat.scatter <- scatter$dat
        attr(res.scatter,"desc") <- dat$desc
        
        meta <- list("dat.scatter"=dat.scatter, "res.scatter"=res.scatter, 
                    "dat.clusters"=dat, "res.clusters"=res)
        
## do auto gating of meta clusters
        G <- res.scatter@K
        childs <- vector("list", G)
        for( g in 1:G ) {
            childs[[g]] <- list("desc"=paste(sep="", "P",g), 
                            "clusters"=.meta.ClustersForScatter(meta,g))
        }
        meta$gating <- list("par"=scatter.subset, "desc"=NULL, "childs"=childs)
        
        par <- rep(TRUE, dat$P)
        par[scatter.subset] <- FALSE
        for( g in 1:G ) {
            meta$gating <- .meta.gating(meta$res.clusters, meta$gating, 
                                    paste(sep="","P",g), which(par), 
                                    c(), iFilter=0)
        }
    }
    else {
        meta <- list("dat.scatter"=NULL, "res.scatter"=NULL,  
                    "dat.clusters"=dat, "res.clusters"=res)
        meta$gating <- list("clusters"=1:res@K, "childs"=c(), 
                    "desc"="all", "partition"=TRUE)

    }
    
    meta
}
## meta.process

##
#   meta.Regression
meta.RegNorm <- function(y, x, method=1, alpha=0.5)
{
    P <- length(y@parameters)
    if( length(x@parameters) != P )
    stop("dimension missmatch")
    
    ys <- y@sigma
    dim(ys) <- c(nrow(y@mu), P*P)
    xs <- x@sigma
    dim(xs) <- c(nrow(x@mu), P*P)
    
    obj <- .C("metaRegNorm",
            as.integer(P), 
            as.integer(y@K), as.double(t(y@mu)), as.double(t(ys)),
            K=as.integer(x@K), M=as.double(t(x@mu)), S=as.double(t(xs)),
            as.integer(method),as.double(alpha),
            package="immunoClust")
    tM <- matrix(obj$M, nrow=obj$K, ncol=P, byrow=TRUE)
    tS <- matrix(obj$S, nrow=obj$K, ncol=(P*P), byrow=TRUE)
    
#colnames(tM) <- colnames(M)
#rownames(tM) <- rownames(M)
    
    list("P"=P, "N"=1, "K"=obj$K, "M"=tM,"S"=tS)

}
##

##
#   meta.Normalize
##
meta.Normalize <- function(P, N, K, W, M, S, G, Z, method=3)
{
    totK <- sum(K)
    
    obj <- .C("metaNormalize", 
            as.integer(P), as.integer(N), as.integer(K),
            W=as.double(W), M=as.double(c(t(M))), S=as.double(c(t(S))), 
            G=as.integer(G), z=as.double(t(Z)), method=as.integer(method), 
            package="immunoCust")
    
    tW <- obj$W
    tM <- matrix(obj$M, nrow=totK, ncol=P, byrow=TRUE)
    tS <- matrix(obj$S, nrow=totK, ncol=(P*P), byrow=TRUE)
    
    colnames(tM) <- colnames(M)
    rownames(tM) <- rownames(M)
    
    list("P"=P, "N"=N, "K"=K, "W"=tW,"M"=tM,"S"=tS)
}
## meta.Normalize


##
#   meta.ME
##
meta.ME <- function(
P, N, K, W, M, S, label, B=100, tol=1e-5, method=20, 
bias=0.25, alpha=0.5, min.class=0
) {
    G <- max(label)
    
    obj <- .Call("immunoC_metaME",
            as.integer(sum(K)), as.integer(P), as.integer(G),
            as.double(W), as.double(c(t(M))), as.double(c(t(S))),
            as.integer(label), as.integer(B), as.double(tol),
            as.integer(method), as.double(bias), as.double(alpha),
            as.integer(min.class))
    
    
    L <- obj$L
# output obj$s to sigma
    sigma <- array(0, c(L, P, P))
    s <- matrix(obj$s, G, P * P, byrow=TRUE)
    for (k in 1:L)
    sigma[k,,] <- matrix(s[k,], P, P, byrow = TRUE)
    
    mu <- matrix(obj$m, G, P, byrow=TRUE)[1:L,]
    dim(mu) <- c(L,P)
    
    z <- matrix(obj$z, sum(K), G, byrow=TRUE)
    z <- as.matrix(z[,1:L])
    
    if( sum(is.na(z)) > 0 ) {
        warning("meta.ME: N/As in Z\n")
        for(row in 1:nrow(z) ) {
            if( sum(is.na(z[row,])) > 0 ) {
                cat("N/As in Z[", row, ",]:", which(is.na(z[row,])), "\n")
            }
        }
        z[is.na(z)] <- 0
    }
#    if( sum( z > 1) > 0 ) {
#        warning("meta.ME: above 1 in Z\n")
#        for(row in 1:nrow(z) ) {
#            if( sum(z[row,]>1) > 0 ) {
#                cat("above 1 in Z[", row, ",]:", which(z[row,]>1), "\n")
#            }
#       }
#       z[z > 1] <- 1
#
#    }
    
# output BIC & ICL
    BIC <- obj$logLike[1]
    ICL <- obj$logLike[2]
    
# output meta-model
    parameters <- colnames(M)
    if( length(parameters) == 0 ) {
        parameters <- paste(sep="", "P", 1:P)
    }
    result <- new("immunoClust", expName="meta.ME", parameters=parameters, 
                    K=L, w=obj$w[1:L], mu=mu, sigma=sigma, 
                    z=z, label=obj$label, 
                    logLike=obj$logLike, BIC=BIC, ICL=ICL,
                    state=obj$history[1:L])
    
    result
    
}
## meta.ME


####
##  meta.Clustering: major iteration loop
###
meta.Clustering <- function(
P, N, K, W, M, S, label=NULL, I.iter=10, B=500, tol=1e-5, 
bias=0.25, alpha=0.5, EM.method=20, 
norm.method=0, norm.blur=2, norm.minG=10
) {

    totK <- sum(K)

    tM <- M
    tS <- S
    
    if( is.null(label) ) {
        label <- rep(1, totK)
        G <- 1
    }
    else {
        G <- max(label)
    }
    for( i in 1:(I.iter) ) {
        if( norm.method > 0 ) {
            if( G >= norm.minG) {
                message("meta.Normalize ", i, " of ", I.iter, " with ", res@K,
                        " clusters (blur=", norm.blur, ")")
                nS <- (1+norm.blur) * tS
                n.res <- meta.ME(P, N, K, W, tM, nS, res@label, B=10, tol=tol, 
                                bias=bias, alpha=alpha, method=10 ) 
                d <- meta.Normalize(P, N, K, W, M, S, 
                                n.res@K, n.res@z, method=norm.method)
                tM <- d$M
                tS <- d$S
            }
        }
        
        if( G < 10 ) 
        label <- meta.SubClustering(P, totK, W, tM, tS, label, tol=tol, 
                                    bias=bias*0.5, alpha=alpha, 
                                    EM.method=EM.method)
        else
        label <- meta.SubClustering(P, totK, W, tM, tS, label, tol=tol, 
                                    bias=bias, alpha=alpha,
                                    EM.method=EM.method)
        
        message("Fit Model ", i, " of ", I.iter, " with ", max(label), 
                " clusters")
        res <- meta.ME(P, N, K, W, tM, tS, label, B=B, tol=tol, 
                        bias=bias, alpha=alpha, method=EM.method)
        
        label <- res@label 
        
        if( res@K == G ) break
        G <- res@K
    }
    attr(res, "bias") <- bias
    res
}


meta.SubClustering <- function(
P, N, W, M, S, label, tol=1e-5, bias=0.25, alpha=1.0, EM.method=20
) {
    
    K <- max(label)
    J <- 8
    cutoff <- 0
    
    res_l <- vector("list", K)
    icl_l <- rep(0, K)
    tst_l <- rep(1, K)
    
    for( k in 1:K ) {
        
        inc <- which(label==k)
        if( length(inc) > 1 ) {
            res <- meta.TestSubCluster(P, length(inc), 
                                    W[inc], M[inc,], S[inc,], 
                                    J=min(length(inc),J), tol=tol, 
                                    bias=bias, alpha=alpha, EM.method=EM.method)
        }
        else {
            res <- NULL
        }
        res_l[[k]] <- res

        if( !is.null(res) && length(res) > 1 ) {
            
            icl <- rep(0, length(res))
            for( l in 1:length(res) )
            icl[l] <- res[[l]]@ICL/res[[l]]@K
            
            icl_l[k] <- max(icl)
            l <- which.max(icl)
            tst_l[k] <- l

#            if( res[[l]]@K == 1 ) {
#                message("cluster ", k, " of ", K, " is OK")
#            }
#            else {
#                message("cluster ", k, " of ", K, ": max. ICL ", 
#                    format(icl_l[k], digits=2), " at ", l, 
#                    " with ", res[[l]]@K, " sub-clusters")
#            }
        }
        else {
            icl_l[k] <- cutoff
            tst_l[k] <- 1
        }
        
    } ## for cluster k
    
    off <- 0
    ins <- vector("list",K)
    sK <- 0
    J <- max(16,2*K)
    
    while( J > sK ) {
        k <- which.max(icl_l)
        
        res <- res_l[[k]]
        icl <- icl_l[k]
        l <- tst_l[k]
        
        if( is.null(res) ) {
            break
        }
        
#if( res[[l]]@K > 1 ) {
#           message("cluster ", k, " has ", res[[l]]@K, " sub-clusters at ", 
#               l, ", ICL=", format(icl, digits=2))
#       }
        
        icl_l[k] <- cutoff
        
        res <- res[[l]]
        
        if( icl <= cutoff )
        break
        
        if( (res@K>1) ) {
            ins[[k]] <- res
        }
        
        sK <- 0
        for(i in 1:K ) 
        if( !is.null(ins[[i]]) ) 
        sK <- sK + (ins[[i]])@K
    }
    
    for( k in 1:K) if( !is.null(ins[[k]]) ) {
#message("split cluster ", k, " into ", (ins[[k]])@K, " sub-clusters")
        label[label==k] <- ins[[k]]@label + max(label)
    }  
    
    label
}
### meta.SubClustering

meta.TestSubCluster<-
function(P, N, W, M, S, J=8, B=500, tol=1e-5, bias=0.5, alpha=1.0, 
        EM.method=2, HC.samples=2000) 
{
#message("meta.TestSubCluster: ", "N=", N, " J=", J, 
#           " dim=", paste(sep="", dim(M), collapes=","))
    if( J > N ) {
        message("\t???", J, "sub >", N, "clusters")
    }
    else {
        
        parameters <- colnames(M)
        if( is.null(parameters) ) {
            parameters <- paste(sep="", "P", 1:P)
        }
        result <- vector("list", J)

        label <- rep(1, N)

        obj <- .Call("immunoC_metaME", 
                as.integer(N), as.integer(P), as.integer(1),
                as.double(W), as.double(c(t(M))), as.double(c(t(S))), 
                label=as.integer(label),
                as.integer(B), as.double(tol), as.integer(EM.method), 
                as.double(bias), as.double(alpha), as.integer(1))

        if( obj$L < 1 ) {
            return(NULL)
        }

        L <- obj$L

# output obj$s to sigma
        sigma <- array(0, c(L, P, P))
        s <- matrix(obj$s, L, P * P, byrow=TRUE)
        for (k in 1:L)
        sigma[k,,] <- matrix(s[k,], P, P, byrow = TRUE)
        
        mu <- matrix(obj$m, L, P, byrow=TRUE)[1:L,]
        dim(mu) <- c(L,P)
        
# output obj$z to Z
        z <- matrix(obj$z, N, 1, byrow=TRUE)
        z <- as.matrix(z[,1:L])
        
# output BIC & ICL
        BIC <- obj$logLike[1]
        ICL <- 0
        logLike <- obj$logLike[3]
        
        
# output
## cat("1: ", obj$logLike, "/n")
        result[[1]] <- new("immunoClust", parameters=parameters, 
                        K=1, N=N, P=P, w=obj$w, 
                        mu=matrix(obj$m, 1, P, byrow=TRUE), 
                        sigma=sigma,label = obj$label,
                        logLike=obj$logLike, BIC=BIC, ICL=ICL)
        
        obj <- NULL
        
# to perform the cluster analysis via EM for each specific number of clusters
        if( J > 1 ) {
##            samples.set = 1:N
## 2016.03.21
## 2016.11.14: till TODO use probability based samples
## probabilities=rowMax( z )
## probes <- 1/(probabilities=rowMax( z ))
            samples.set <- which( apply(S, 1, function(s)
                                    { 
                                    dim(s) <- c(P,P)
                                    d <- det(s)
                                    if( d <= 0 ) return (-300)
                                        log(d)
                                    }) > -100)
            if( length(samples.set) < N ) {
                cat(length(samples.set), "of", N, "clusters used\n")
            }
            
            if( J > length(samples.set) ) {
                return (NULL)
            }
            if( length(samples.set) > HC.samples ) {
                samples.set = sample(samples.set, HC.samples)
                hcPairs <- meta.hclust(P, HC.samples, W[samples.set], 
                                        M[samples.set,], S[samples.set,])
            }
            else {
                HC.samples <- length(samples.set)
                hcPairs <- meta.hclust(P, HC.samples, W[samples.set], 
                                        M[samples.set,], S[samples.set,])
            }
            
            for (k in 2:J) {
                
                label <- rep(0, N)
                
                label[samples.set] <- .clust.hclass(hcPairs, k)
                
                obj <- .Call("immunoC_metaME", 
                        as.integer(N), as.integer(P), as.integer(k),
                        as.double(W), as.double(c(t(M))), as.double(c(t(S))), 
                        label=as.integer(label), 
                        as.integer(B), as.double(tol), as.integer(EM.method), 
                        as.double(bias), as.double(alpha), as.integer(1))

                L <- obj$L

# output obj$s to sigma
                sigma <- array(0, c(L, P, P))
                s <- matrix(obj$s, k, P * P, byrow=TRUE)
                if( L > 0 ) {
                    for (l in 1:L)
                        sigma[l,,] <- matrix(s[l,], P, P, byrow = TRUE)

                    mu <- matrix(obj$m, k, P, byrow=TRUE)[1:L,]
                    dim(mu) <- c(L,P)
                }   
                else {
                    mu <- matrix(0, nrow=0, ncol=P)
                }
# output BIC & ICL
                BIC <- obj$logLike[1]
                ICL <- obj$logLike[3] - logLike
                
# outp
                result[[k]] <- new("immunoClust", parameters=parameters, 
                                K=L, N=N, P=P, w=obj$w, mu=mu, sigma=sigma,
                                label = obj$label,
                                logLike=obj$logLike, BIC=BIC, ICL=ICL)
                
                obj <- NULL
                
                
            } ## for k
        } ## if J>1
        
        result
        
    }
    
}
### meta.TestSubCluster

###
# meta.hclust
###
meta.hclust <- function(P, N, W, M, S) 
{
    partition <- 1:N
    attr(partition, "unique") <- N
    
    obj <- .C("metaHC",
            li=integer(N-1), lj=integer(N-1), crit=double(N-1),
            as.integer(N), as.integer(P), 
            as.double(W), c(t(M)), c(t(S)),
            package="immunoClust")
    merge <- .clust.hpairs2(obj$li, obj$lj)
    structure(t(merge), initialPatition=partition, change=obj$crit,
                dimensions=c(N,P), modelName="meta",
                call = match.call())
}
### meta.hclust
