####
# meta.process
####
meta.process <- function(
exp, dat.subset=c(), meta.iter=10, tol=1e-5, meta.bias=0.2, 
meta.alpha=.5, norm.method=0, norm.blur=2, norm.minG=10
) {
    dat <- meta.exprs(exp, sub=dat.subset)
    
    res <- meta.Clustering(dat$P, dat$N, dat$K, 
                        dat$clsEvents, dat$M, dat$S, 
                        bias=meta.bias, sub.thres=0, I.iter=meta.iter, B=50, 
                        tol=tol, EM.method=20, alpha=meta.alpha,
                        norm.method=norm.method, norm.blur=norm.blur, 
                        norm.minG=norm.minG)

    dat.norm <- dat
    if( norm.method > 0 ) {
        dat.norm <- meta.Normalize(dat$P, dat$N, dat$K, dat$clsEvents, 
                                dat$M, dat$S, res@K, res@z, 
                                method=norm.method)
    }
    
    attr(res, "desc") <- dat$desc
    if( is.null(dat.subset) )
        dat.subset <- seq_len(npar(res))
    attr(res, "trans.a") <- apply( vapply(exp,function(x) x@trans.a[dat.subset],
                                rep(0.01,npar(res))), 1, mean)
    attr(res, "trans.b") <- apply( vapply(exp,function(x) x@trans.b[dat.subset],
                                rep(0.0,npar(res))), 1, mean)
    attr(res, "limits") <- attr(exp[[1]], "limits")[,dat.subset]

    meta <- list("dat.scatter"=NULL, "res.scatter"=NULL,  
                    "dat.clusters"=dat, "res.clusters"=res)
    meta$gating <- list("clusters"=seq_len(res@K), "childs"=c(),
                    "desc"="all", "partition"=TRUE)

    
    meta$gating$pscales <- Default_Scales(attr(res, "trans.a"),
                                        attr(res, "limits"))
    
    class(meta) <- "immunoMeta"
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
    
    if( method == 23 ) {
        ## 2018.12.06:
        ## method==23: K[1] cluster are not re-labeled
        ## min.class parameter is re-used for this method
        min.class <- K[1]
    }
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
    for (k in seq_len(L))
    sigma[k,,] <- matrix(s[k,], P, P, byrow = TRUE)
    
    mu <- matrix(obj$m, G, P, byrow=TRUE)[seq_len(L),]
    dim(mu) <- c(L,P)
    
    z <- matrix(obj$z, sum(K), G, byrow=TRUE)
    z <- as.matrix(z[,seq_len(L)])
    
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
        parameters <- paste(sep="", "P", seq_len(P))
    }
    result <- new("immunoClust", expName="meta.ME", parameters=parameters, 
                    K=L, P=P, w=obj$w[seq_len(L)], mu=mu, sigma=sigma, 
                    z=z, label=obj$label, 
                    logLike=obj$logLike, BIC=BIC, ICL=ICL,
                    state=obj$history[seq_len(L)])
    
    result
    
}
## meta.ME


####
##  meta.Clustering: major iteration loop
###
meta.Clustering <- function(
P, N, K, W, M, S, label=NULL, I.iter=10, B=500, tol=1e-5, 
bias=0.25, sub.thres=bias, alpha=0.5, EM.method=20, HC.samples=2000,
norm.method=0, norm.blur=2, norm.minG=10, verbose=FALSE
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
    
    ## 2022.10.18:
    res <- new("immunoClust", expName="metaClust", 
            K=G, N=totK, P=P, label=label)
    
    subEM.method <- EM.method
    ### 2018.12.06: trial EM.method=23 to fix cluster of first sample
    ### sub-clustering is performed with EM.method=20
    if( EM.method == 23 )
        subEM.method <- 20
    for( i in seq_len(I.iter) ) {
        if( norm.method > 0 ) {
            if( G >= norm.minG) {
                #message("meta.Normalize ", i, " of ", I.iter, " with ", res@K,
                #        " clusters (blur=", norm.blur, ")")
                nS <- (1+norm.blur) * tS
                n.res <- meta.ME(P, N, K, W, tM, nS, res@label, B=10, tol=tol, 
                                bias=bias, alpha=alpha, method=10 ) 
                d <- meta.Normalize(P, N, K, W, M, S, 
                                n.res@K, n.res@z, method=norm.method)
                tM <- d$M
                tS <- d$S
            }
        }
        
        if( res@K < 10 )
        sub_bias <- bias*0.5
        else
        sub_bias <- bias
        label <- meta.SubClustering(res, P, totK, W, tM, tS,
                                    tol=tol,
                                    bias=sub_bias, thres=sub.thres, alpha=alpha,
                                    EM.method=subEM.method, HC.samples=HC.samples,
                                    verbose=verbose)
        
        message("Fit Model ", i, " of ", I.iter, " with ",
                length(unique(label)), " clusters")
        res <- meta.ME(P, N, K, W, tM, tS, label, B=B, tol=tol, 
                        bias=bias, alpha=alpha, method=EM.method)
        message("=> results in  ", res@K, " clusters")
        label <- res@label 
        
        if( res@K == G ) break
        G <- res@K
    }
    attr(res, "bias") <- bias
    res
}


#meta.SubClustering <- function(
#P, N, W, M, S, label, tol=1e-5, bias=0.25, alpha=1.0, EM.method=20
meta.SubClustering <- function(
x, P, N, W, M, S, tol=1e-5, bias=0.25, thres=bias, alpha=1.0, 
EM.method=20, HC.samples=2000,
verbose=FALSE
) {
    
    label <- x@label
    K <- x@K
    #K <- max(label)
    
    J <- 8
  
    
    ## 2022.010.10: find a balance to introduce more clusters
    ##thres <- bias
    icl_thres <- (P*(P+1)/2 + P)*log(sum(W))*0.5*thres
    #message("ICL thres ",  format(icl_thres, digits=2) )
    #cutoff <- min(0, icl_thres)
    cutoff <- 0
    
    res_l <- vector("list", K)
    icl_l <- rep(0, K)
    tst_l <- rep(1, K)
    
    for( k in seq_len(K) ) {
        
        inc <- which(label==k)
        
        if( length(inc) > 1 ) {
            res <- meta.TestSubCluster(x,
                    P, length(inc), W[inc], M[inc,], S[inc,],
                    J=min(length(inc),J), tol=tol,
                    bias=bias, alpha=alpha, 
                    EM.method=EM.method, HC.samples=HC.samples )
                    
            if( is.null(res) && verbose ) {
                message("cluster ", k, " null TestSubCluster")
            }
        }
        else {
            res <- NULL
        }
        
        
        res_l[[k]] <- res

        if( !is.null(res) && length(res) > 1 ) {
            
            icl <- rep(0, length(res))
            for( l in seq_len(length(res)) )
            icl[l] <- res[[l]]@ICL/res[[l]]@K
            
            icl[1] <- cutoff
            
            icl_l[k] <- max(icl)
            l <- which.max(icl)
            tst_l[k] <- l
            
            
            #if( verbose ) {
            #    ccl <- rep(0, length(res))
            #    for( l in seq_len(length(res)) )
            #    ccl[l] <- res[[l]]@logLike[3]
            #
            #    message("TestSubCluster ", bias, ": ", k, " ", length(inc), ">>",
            #        length(res), ">>", paste(round(icl), collapse=","), "<<", paste(round(ccl), collapse=",") )
            #}
           
           
        }
        else {
            icl_l[k] <- cutoff
            tst_l[k] <- 1
        }
        
        #if( verbose ) {
        #    message("cluster", k, ": icl", icl_l[k], "(", tst_l[k], ")")
        #}
        
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
            if( verbose ) {
                message(J, "/", sK, ": break with null at", k)
            }
            break
        }
        
        if( res[[l]]@K > 1 && verbose ) {
            message(J, "/", sK, ": cluster ", k, " has ", res[[l]]@K, 
                " sub-clusters at ", l, ", ICL=", format(icl, digits=2), "<>",
                format(icl_thres, digits=2))
        }
        
        icl_l[k] <- cutoff
        
        res <- res[[l]]
        
        if( icl <= cutoff ) {
            if( verbose ) {
                message("break at ", k, " ", icl, " ", cutoff)
            }
            break
        }
        
        ## 2022.10.17: respect cluster costs
        if( (res@K>1) && (icl>icl_thres) ) {
            ins[[k]] <- res
        }
        
        sK <- 0
        for(i in seq_len(K) )
        if( !is.null(ins[[i]]) ) 
        sK <- sK + (ins[[i]])@K
    } # while J > sK
    
    for( k in seq_len(K) ) if( !is.null(ins[[k]]) ) {
        #message("split cluster ", k, " into ", (ins[[k]])@K, " sub-clusters")
## 2018.12.05, try do better re-labeling with the aim:
## first label==k cluster and all of them with same sub-cluster label
## keep label k, the rest get's a new label
## should have the consequence, that the cluster label of the first experiment
## is not changed if they all distinct??? (required for EM.method=23)
        inc <- which(label==k)
        lnc <- (ins[[k]]@label==ins[[k]]@label[1])
        knc <- inc[lnc]
        label[label==k] <- ins[[k]]@label + max(label)
        label[knc] <- k
    }  #for cluster k
    
    label
}
### meta.SubClustering

meta.TestSubCluster<-
function(x, P, N, W, M, S, J=8, B=500, tol=1e-5, bias=0.5, alpha=1.0,
        EM.method=2, HC.samples=2000)
{
#message("meta.TestSubCluster: ", "N=", N, " J=", J, 
#           " dim=", paste(sep="", dim(M), collapes=","))
    if( J > N ) {
        warning("\t???", J, "sub >", N, "clusters")
        return (NULL)
    }

    K <- x@K
    sumW <- sum(W)
    
    parameters <- colnames(M)
    if( is.null(parameters) ) {
        parameters <- paste(sep="", "P", seq_len(P))
    }
    result <- vector("list")
    
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
    for (k in seq_len(L))
    sigma[k,,] <- matrix(s[k,], P, P, byrow = TRUE)
        
    mu <- matrix(obj$m, L, P, byrow=TRUE)[seq_len(L),]
    dim(mu) <- c(L,P)
        
# output obj$z to Z
    z <- matrix(obj$z, N, 1, byrow=TRUE)
    z <- as.matrix(z[,seq_len(L)])
        
# output BIC & ICL
    BIC <- obj$logLike[1]
    ICL <- 0
    logLike <- obj$logLike[3]
        
        
# output
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
        
        
        if( J > length(samples.set) ) {
            ##message( J, "<>", samples.set, " N=", N)
            return (result)
        }
        
        if( length(samples.set) > HC.samples ) {
            samples.set <- sample(samples.set, HC.samples)
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
                for (l in seq_len(L))
                    sigma[l,,] <- matrix(s[l,], P, P, byrow = TRUE)
                    
                mu <- matrix(obj$m, k, P, byrow=TRUE)[seq_len(L),]
                dim(mu) <- c(L,P)
            }
            else {
                mu <- matrix(0, nrow=0, ncol=P)
            }
# output BIC & ICL
            BIC <- obj$logLike[1]
                ## 2022.10.10: .icl_delta ???
            #ICL <- obj$logLike[3] - logLike
            ICL <- obj$logLike[3] - logLike - .icl_delta(sumW, P, K, L)*bias
                
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
### meta.TestSubCluster

###
# meta.hclust
###
meta.hclust <- function(P, N, W, M, S) 
{
    partition <- seq_len(N)
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
