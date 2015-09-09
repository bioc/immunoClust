####
# meta.process
####
meta.process <- function(
exp, dat.subset=c(), meta.iter=10, meta.bias=0.2, 
meta.alpha=.5, meta.normalize=FALSE, norm.degree=1, 
scatter.subset=c(1,2), scatter.bias=0.25,scatter.prior=6
) {
    dat <- meta.exprs(exp, sub=dat.subset)
    
    res <- meta.Clustering(dat$P, dat$N, dat$K, dat$clsEvents, dat$M, dat$S, 
                        bias=meta.bias, I.iter=meta.iter, EM.method=20, 
                        alpha=meta.alpha, normalize=meta.normalize, 
                        norm.degree=norm.degree)

    dat.norm <- dat
    if( meta.normalize ) {
        dat.norm <- meta.Normalize(dat$P, dat$N, dat$K, dat$clsEvents, 
                                dat$M, dat$S, res@K, res@z, degree=norm.degree)
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
        
        childs <- vector("list", 1)
        childs[[1]] <- list("desc"="P1", "clusters"=1:(res@K))
        meta$gating <- list("par"=scatter.subset, "desc"=NULL, "childs"=childs)
        
        meta$gating <- .meta.gating(meta$res.clusters, meta$gating, 
                                "P1", 1:(dat$P), c(), iFilter=0)
    }
    
    meta
}
## meta.process

##
#   meta.Scale
##
meta.Scale <- function(P, N, K, W, M, S, method=5)
{
    totK <- sum(K)
    label <- rep(1, totK)
    obj <- .C("metaScale", 
            as.integer(P), as.integer(N), as.integer(K),
            W=as.double(W), M=as.double(c(t(M))), S=as.double(c(t(S))), 
            label=as.integer(label), as.integer(method), 
            package="immunoClust")

    tW <- obj$W
    tM <- matrix(obj$M, nrow=totK, ncol=P, byrow=TRUE)
    tS <- matrix(obj$S, nrow=totK, ncol=(P*P), byrow=TRUE)

    colnames(tM) <- colnames(M)
    rownames(tM) <- rownames(M)

    list("P"=P, "N"=N, "K"=K, "W"=tW,"M"=tM,"S"=tS)
}
## meta.Scale

##
#   meta.Mormalize
##
meta.GPA <- function(P, N, K, W, M, S, G, Z)
{
    totK <- sum(K)
    groups <- rep(1, totK)
    
    obj <- .C("metaGPA", 
            as.integer(P), as.integer(N), as.integer(K),
            W=as.double(W), M=as.double(c(t(M))), S=as.double(c(t(S))), 
            G=as.integer(G), z=as.double(t(Z)), groups=as.integer(groups), 
            package="immunoCust")
    
    tW <- obj$W
    tM <- matrix(obj$M, nrow=totK, ncol=P, byrow=TRUE)
    tS <- matrix(obj$S, nrow=totK, ncol=(P*P), byrow=TRUE)
    
    colnames(tM) <- colnames(M)
    rownames(tM) <- rownames(M)
    
    list("P"=P, "N"=N, "K"=K, "W"=tW,"M"=tM,"S"=tS)
}
meta.Normalize <- function(P, N, K, W, M, S, G, Z, degree=1)
{
    totK <- sum(K)
    
    obj <- .C("metaNormalize", 
            as.integer(P), as.integer(N), as.integer(K),
            W=as.double(W), M=as.double(c(t(M))), S=as.double(c(t(S))), 
            G=as.integer(G), z=as.double(t(Z)), degree=as.integer(degree), 
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
P, N, K, W, M, S, label, B=500, tol=1e-5, method=20, 
bias=0.25, alpha=0.5, min.class=0
) {
    G <- max(label)
#    obj <- .C("metaME", 
#            as.integer(P), as.integer(N), as.integer(K),
#            as.double(W), as.double(c(t(M))), as.double(c(t(S))), 
#            G=as.integer(G), w=double(G), m=double(G*P), s=double(G*P*P),
#            label=as.integer(label), logLike=double(3), 
#            history=as.integer(1:G),
#            as.integer(B), as.double(tol), 
#            as.integer(method), as.double(bias), 
#            as.double(alpha), as.integer(min.class),
#            package="immunoClust")
    
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
        z[is.na(z)] <- 0
    }
    
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
P, N, K, W, M, S, I.iter=10, B=500, tol=1e-5, 
bias=0.25, alpha=0.5, EM.method=20, 
normalize=FALSE, norm.degree=1, norm.minG=10
) {
    
    totK <- sum(K)
    tM <- M
    tS <- S
    if( normalize ) {
        d <- meta.Scale(P, N, K, W, M, S, method=5)
        tM <- d$M
        tS <- d$S
    }
    
    label <- rep(1, totK)
    G <- 1
    for( i in 1:(I.iter) ) {
        if( normalize && G >= norm.minG) {
            d <- meta.Normalize(P, N, K, W, M, S, 
                                res@K, res@z, degree=norm.degree)
            tM <- d$M
            tS <- d$S
        }
        
        
        if( G < 10 ) 
        label <- meta.SubClustering(P, totK, W, tM, tS, label, tol=tol, 
                                    bias=bias*0.5, alpha=alpha, 
                                    EM.method=EM.method)
        else
        label <- meta.SubClustering(P, totK, W, tM, tS, label, tol=tol, 
                                    bias=bias, alpha=alpha,
                                    EM.method=EM.method)
        
        message("meta.Clustering ", i, " of ", I.iter)
        res <- meta.ME(P, N, K, W, tM, tS, label, B=B, tol=tol, 
                        bias=bias, alpha=alpha, method=EM.method)
        
        label <- res@label 
        
        if( res@K == G ) break
        G <- res@K
    }
    
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
            message("cluster ", k, " of ", K, ": TestSubCluster")
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
            
            if( res[[l]]@K == 1 ) {
                message("cluster ", k, " of ", K, " is OK")
            }
            else {
                message("cluster ", k, " of ", K, ": max. ICL ", 
                    format(icl_l[k], digits=2), " at ", l, 
                    " with ", res[[l]]@K, " sub-clusters")
            }
        }
        else {
            icl_l[k] <- cutoff
            tst_l[k] <- 1
        }
        
    } ## for cluster k
    
    off <- 0
    ins <- vector("list",K)
    sK <- 0
    J <- max(8,2*K)
    
    while( J > sK ) {
        k <- which.max(icl_l)
        
        res <- res_l[[k]]
        icl <- icl_l[k]
        l <- tst_l[k]
        
        if( is.null(res) ) {
            break
        }
        
        if( res[[l]]@K > 1 ) {
            message("cluster ", k, " has ", res[[l]]@K, " sub-clusters at ", 
                l, ", ICL=", format(icl, digits=2))
        }
        
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
        message("split cluster ", k, " into ", (ins[[k]])@K, " sub-clusters")
        label[label==k] <- ins[[k]]@label + max(label)
    }  
    
    label
}
### meta.SubClustering

meta.TestSubCluster<-
function(P, N, W, M, S, J=8, B=500, tol=1e-5, bias=0.5, alpha=1.0, 
        EM.method=2, HC.samples=2000) 
{
    message("meta.TestSubCluster: ", "N=", N, " J=", J, 
            " dim=", paste(sep="", dim(M), collapes=","))
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

#        obj <- .C("metaME", 
#                as.integer(P), as.integer(1), as.integer(N),
#                as.double(W), as.double(c(t(M))), as.double(c(t(S))), 
#                G=as.integer(1), w=double(1), m=double(1*P), s=double(1*P*P),
#                label=as.integer(label), logLike=double(3), 
#                history=as.integer(1),
#                as.integer(B), as.double(tol), as.integer(EM.method), 
#                as.double(bias), as.double(alpha), as.integer(0),
#                package="immunoClust")
#        
#        if( obj$G < 1 ) return(NULL)
#        
#        L <- obj$G

        obj <- .Call("immunoC_metaME", 
                as.integer(N), as.integer(P), as.integer(1),
                as.double(W), as.double(c(t(M))), as.double(c(t(S))), 
                label=as.integer(label),
                as.integer(B), as.double(tol), as.integer(EM.method), 
                as.double(bias), as.double(alpha), as.integer(0))

        if( obj$L < 1 ) return(NULL)

        L <- obj$L

# output obj$s to sigma
        sigma <- array(0, c(L, P, P))
        s <- matrix(obj$s, L, P * P, byrow=TRUE)
        for (k in 1:L)
        sigma[k,,] <- matrix(s[k,], P, P, byrow = TRUE)
        
        mu <- matrix(obj$m, L, P, byrow=TRUE)[1:L,]
        dim(mu) <- c(L,P)
        
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
            samples.set = 1:N
            if( N > HC.samples ) {
                samples.set = sample(1:N, HC.samples)
                hcPairs <- meta.hclust(P, HC.samples, W[samples.set], 
                                        M[samples.set,], S[samples.set,])
            }
            else {
                hcPairs <- meta.hclust(P, N, W, M, S)
            }
            
            for (k in 2:J) {
                
                label <- rep(0, N)
                
                label[samples.set] <- .clust.hclass(hcPairs, k)
                
#                obj <- .C("metaME", 
#                        as.integer(P), as.integer(1), as.integer(N),
#                        as.double(W), as.double(c(t(M))), as.double(c(t(S))), 
#                        G=as.integer(k), w=double(k), 
#                        m=double(k*P), s=double(k*P*P),
#                        label=as.integer(label), logLike=double(3), 
#                        history=1:k,
#                        as.integer(B), as.double(tol), as.integer(EM.method), 
#                        as.double(bias), as.double(alpha), as.integer(0),
#                        package="immunoClust")
#                                
#                L <- obj$G

                obj <- .Call("immunoC_metaME", 
                        as.integer(N), as.integer(P), as.integer(k),
                        as.double(W), as.double(c(t(M))), as.double(c(t(S))), 
                        label=as.integer(label), 
                        as.integer(B), as.double(tol), as.integer(EM.method), 
                        as.double(bias), as.double(alpha), as.integer(0))

                L <- obj$L

# output obj$s to sigma
                sigma <- array(0, c(L, P, P))
                s <- matrix(obj$s, k, P * P, byrow=TRUE)
                for (l in 1:L)
                sigma[l,,] <- matrix(s[l,], P, P, byrow = TRUE)
                
                mu <- matrix(obj$m, k, P, byrow=TRUE)[1:L,]
                dim(mu) <- c(L,P)
                
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
