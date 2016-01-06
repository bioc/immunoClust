##
#    meta.scatter.gating
##
.meta.scatter.prior.6 <- function(dat)
{
    M <- dat$M
    a <- .separation.maximization(M[,1], more=TRUE)
    b <- .separation.maximization(M[,2], more=TRUE)
    N <- nrow(M)
    label <- rep(0, N)
    
# bkg
    label[M[,1] < a[1] & M[,2] < b[1]] <- 6
# ly
    label[M[,1] > a[1] & M[,1] < a[3] & M[,2] < b[1] ] <- 1
# gr
    label[M[,1] > a[1] & M[,1] < a[3] & M[,2] > b[2] ] <- 3
# mo
    label[M[,1] > a[2] & M[,1] < a[3] & M[,2] > b[1] & M[,2] < b[3] ] <- 2
    
# gr dupl
    label[M[,1] > a[3] & M[,2] > b[3] ] <- 4
# mo dupl
    label[M[,1] > a[3] & M[,2] > b[2] & M[,2] < b[3] ] <- 5
# ly dupl
#    label[M[,1] > a[3] & M[,2] < b[2]  ] <- 6
    
    gates <- array(NA, 2*3)
    dim(gates) <- c(1, 2, 3)
    gates[1,1,1:length(a)] <- a
    gates[1,2,1:length(b)] <- b
    
    
    list("label"=label, "names"=paste(sep="", "P", 1:6), gates=gates)
}
.meta.scatter.prior.7 <- function(dat)
{
    M <- dat$M
    a <- .separation.maximization(M[,1], more=TRUE)
    b <- .separation.maximization(M[,2], more=TRUE)
    N <- nrow(M)
    label <- rep(0, N)
    
# bkg
    label[M[,1] < a[1] & M[,2] < b[1]] <- 7
# ly
    label[M[,1] > a[1] & M[,1] < a[3] & M[,2] < b[1] ] <- 1
# gr
    label[M[,1] > a[1] & M[,1] < a[3] & M[,2] > b[2] ] <- 3
# mo
    label[M[,1] > a[2] & M[,1] < a[3] & M[,2] > b[1] & M[,2] < b[3] ] <- 2
    
# gr dupl
    label[M[,1] > a[3] & M[,2] > b[3] ] <- 4
# mo dupl
    label[M[,1] > a[3] & M[,2] > b[2] & M[,2] < b[3] ] <- 5
# ly dupl
    label[M[,1] > a[3] & M[,2] < b[2]  ] <- 6
    
    gates <- array(NA, 2*3)
    dim(gates) <- c(1, 2, 3)
    gates[1,1,1:length(a)] <- a
    gates[1,2,1:length(b)] <- b
    
    
    list("label"=label, "names"=paste(sep="", "P", 1:7), gates=gates)
}

.meta.scatter.gating <- function(
exp, sub=c(1,2), more=c(TRUE,TRUE), scatter.prior=6, 
EM.method=1, EM.bias=0.25, EM.alpha=.5
) {
    
    dat <- meta.exprs(exp, sub=sub)
    
    if( scatter.prior == 7 ) {
        meta <- .meta.scatter.prior.7(dat)
    }
    else {
        meta <- .meta.scatter.prior.6(dat)
    }
    cat("Scatter Gating\n")
    try(res <- meta.ME(dat$P, dat$N, dat$K, dat$clsEvents, dat$M, dat$S, 
                    meta$label, method=EM.method, bias=EM.bias, alpha=EM.alpha))
    
    meta$label <- res@label
    G <- res@K
    gates <- array(NA, G*length(sub)*3)
    dim(gates) <- c(G, length(sub),3)
    childs <- vector("list", G)
    for( g in 1:G ) {
        gates[g,,] <- as.matrix(meta$gates[1,,])
        childs[[g]] <- list("desc"=paste(sep="", "P",g), "clusters"=g)
        
    }
    gating <- list("par"=sub, "desc"=NULL, 
                "gates"=meta$gates[1,,], "childs"=childs)

    meta$gates <- gates
    meta$par <- sub
    meta$scatter <- res
    list("dat"=dat, "res"=res, "meta"=meta, "gates"=gates, "gating"=gating)
}
## meta.scatter.gating


.meta.getGate <- function(gating, pattern)
{
    res <- NULL
    N <- length(pattern)
    
    if( N == 0 ) {
        res <- gating
    }
    else 
    if( !is.null(gating$childs) ) {
        for( i in 1:length(gating$childs) ) {
            if( gating$childs[[i]]$desc == pattern[1] ) {
                if( N == 1 ) {
                    res <- gating$childs[[i]]
                }
                else {
                    res <- .meta.getGate(gating$childs[[i]], pattern[2:N])
                }
            }
        }
    }
    
    res
    
}

.meta.setGate <- function(gating, pattern, gate)
{
    N <- length(pattern)
    
    if( N == 0 ) {
        gating <- gate
    }
    else 
    if( !is.null(gating$childs) ) {
        for( i in 1:length(gating$childs) ) {
            if( gating$childs[[i]]$desc == pattern[1] ) {
                if( N == 1 ) {
                    gating$childs[[i]] <- gate
                }
                else {
                    gating$childs[[i]] <- 
                        .meta.setGate(gating$childs[[i]], pattern[2:N], gate)
                }
            }
        }
    }
    
    gating
    
}


##
# meta.gating
##
.meta.gating <- 
function(res, gating, pattern, par, par.order=NULL, 
        iFilter=0, aFac=5/3, bFac=3/2)
{
    gate <- .meta.getGate(gating,pattern)
    
    inc <- rep(FALSE, res@K)
    inc[gate$clusters] <- TRUE
    
    post <- .meta.posterior.gating(res, inc, par, par.order=par.order, 
                                iFilter=iFilter, aFac=aFac, bFac=bFac)
    post$desc <- gate$desc
    post$clusters <- gate$clusters
    
    gating <- .meta.setGate(gating, pattern, post)
    
    gating
}
## meta.gating

##
# meta.posterior.gating
##
.meta.posterior.gating <- 
function(res, inc, par, par.order=c(), iFilter=0, aFac=7/6, bFac=3/2, pop="")
{
    desc <- res@desc
    
    ret <- NULL
    if( sum(inc) > 0  ) {
        if( length(par) > 0 ) {
            gates <- array(NA, 3*length(par))
            dim(gates) <- c(length(par),3)
            score <- rep(0, length(par))
            
## filter meta-clusters with only one experiment clusters
            unc <- inc
            for( i in 1:length(inc) ) {
                if( inc[i] ) {
                    if( sum(res@label==i) <= iFilter ) {
                        unc[i] <- FALSE
                    }
                }
            }
            
            if( sum(unc) > 0 ) {
                for( i in 1:(length(par)) ) {
                    p <- par[i]
                    M <- res@mu[unc, p]
                    S <- res@sigma[unc, p, p]
                    g <- .classification.maximization(M, S, 
                                                    aFac=aFac, bFac=bFac)
                    if( is.na(g[1]) ) {
                        score[i] <- bFac
                    }
                    else {
                        score[i] <- .classification.score(M, S, g)
                        if( g[1] == g[2] ) {
## less gates has better score
                            score[i] <- score[i]*aFac
                        }
                        lbl <- rep(2, sum(unc))
                        lbl[ M <= g[1] ] <- 1
                        lbl[ M >= g[2] ] <- 3
                        sm <- .classification.separation(res@mu[unc,], 
                                                    res@sigma[unc,,], lbl) 
                    }
                    gates[i,] <- c(g,NA)
                }
            }
            
            i <- which.max(score)

            nogate <- is.na(gates[i,1])
            
            childs <- vector("list")
####
            if( !is.na(gates[i,1]) ) {
                
                zsuf <- c("neg", "dim", "pos")
                
                if( gates[i,1] != gates[i,2] ) {
                    p <- par[i]
                    g <- gates[i,]
                    
                    par <- par[ par!= p ]
                    
                    c <- 0
                    M <- res@mu
                    
## neg
                    lbl <- inc & M[,p] <= g[1]
                    child <- .meta.posterior.gating(res, lbl, par, 
                                par.order=par.order, iFilter=iFilter, 
                                aFac=aFac, bFac=bFac, 
                                pop=paste(sep="", pop, "_", desc[p], "neg"))
                    if( !is.null(child) ) {
                        child$desc <- paste(sep="", desc[p], "neg")
                        c <- c+1
                        childs[[c]] <- child
                    }
                    
## pos
                    lbl <- inc & (M[,p] >= g[2])
                    child <- .meta.posterior.gating(res, lbl, par, 
                                par.order=par.order, iFilter=iFilter, 
                                aFac=aFac, bFac=bFac, 
                                pop=paste(sep="", pop, "_", desc[p], "pos"))
                    if( !is.null(child) ) {
                        child$desc <- paste(sep="", desc[p], "pos")
                        c <- c+1
                        childs[[c]] <- child
                    }
## dim
                    lbl <- inc & (M[,p] > g[1]) & (M[,p] < g[2])
                    child <- .meta.posterior.gating(res, lbl, par, 
                                par.order=par.order, iFilter=iFilter, 
                                aFac=aFac, bFac=bFac, 
                                pop=paste(sep="", pop, "_", desc[p], "dim"))
                    if( !is.null(child) ) {
                        child$desc <- paste(sep="", desc[p], "dim")
                        c <- c+1
                        childs[[c]] <- child
                    }
                                        
                }
                else 
                if( score[i] < 3 ) {                    
                    nogate <- TRUE
                }
                else {
                    
### sort by score
                    map <- match(length(score):1, rank(score, ties.method="f"))
                    
                    par <- par[map]
                    gates[1:length(map),] <- gates[map,]
                    score <- score[map]
                    i <- which.max(score)
                    
###                
                    ipar <- !is.na(gates[,1]) & (gates[,1]==gates[,2]) & 
                            (score >= max(3.0,score[i]/aFac))
                    
                    zpar <- which( ipar )
                    npar <- which( !ipar )
                    
## par.order
                    if( !is.null(par.order) ) {
                        opar <- c()
                        for( p in par.order ) {
                            i <- match(p, par[zpar])
                            if( !is.na(i) ) {
                                opar <- c(opar, zpar[i])
                            }
                        }
                        for( p in zpar ) {
                            if( !match(p, opar) ) {
                                opar <- c(opar, p)
                            }
                        }
                        
                        zpar <- opar
                    }
                    
                    c <- 0
                    M <- res@mu
                    
                    build.childs <- function(zbl, jpar, pop="")
                    {
                        
                        childs <- vector("list")
                        c <- 0
                        p <- par[jpar]            ## index in par
                        g <- gates[jpar,]
                        
                        for( z in 1:3 ) {

                            lbl <- zbl
#                
                            if( z == 1 ) {
## neg
                                lbl <- zbl & (M[,p] <= g[1])
                            }
                            else
                            if( z == 3 ) {
## pos
                                lbl <- zbl & (M[,p] >= g[2])
                            }
                            else { 
## dim
                                lbl <- zbl & (M[,p] > g[1]) & (M[,p] < g[2])
                            }
                            
                            znam <- paste(sep="", desc[p], zsuf[ z ])
                            
                            if( sum(lbl) > 0) {
                                child <- NULL
                                if( jpar < length(zpar) ) {
                                    pchlds <- build.childs(lbl, jpar+1, 
                                        pop=paste(sep="", pop, "_", znam))
                                    child <- list("par"=par[jpar+1], 
                                                "desc"=desc[p], 
                                                "gates"=gates[jpar+1,], 
                                                "childs"=pchlds, 
                                                "clusters"=which(lbl))
                                }
                                else {
                                    child <- .meta.posterior.gating(
                                        res, lbl, par[npar], 
                                        par.order=par.order, iFilter=iFilter, 
                                        aFac=aFac, bFac=bFac, 
                                        pop=paste(sep="", pop, "_", znam))
                                    
                                }
                                child$desc <- znam
                                c <- c+1
                                childs[[c]] <- child
                            }
                        }
                        
                        childs
                    }
                    
                    childs <- build.childs(inc, 1, pop=pop)
                    p <- par[ 1 ]
                    g <- gates[1, ]
                    
                }
            }
####
            if( nogate ) {
## leave node 
                ret <- list("desc"="", "clusters"=which(inc))
            }
            else {
## branch node            
                ret <- list("par"=p, "desc"=desc[p], "gates"=g, 
                        "childs"=childs, "clusters"=which(inc))
            }
        }
        else {
## leave node: length(par)==0
            ret <- list("desc"="", "clusters"=which(inc))
        }
## length(par) > 0 
        
    }
## sum(inc) > 0
    
## cat("gating<<<<<", pop, "\n")
    
    ret
}

###

#meta.expForClusters <- function(meta, clusters)
#{
#    N <- meta$dat.clusters$N
#    K <- meta$dat.clusters$K
#
#    res <- rep(0,N)
#    for( a in clusters ) {
#        label <- which(meta$res.clusters@label == a)
#    
#        for( l in label ) {
#            k <- 0
#            for( i in 1:N ) {
#                k <- k + K[i]
#                if( l <= k ) {
#                    res[i] <- res[i] + 1
#                    break
#                }
#            }
#        }
#    }
#    res
#}


.meta.ClustersForScatter <- function(meta, scatter, filter=0)
{
    ret <- c()
    scatter_label <- meta$res.scatter@label
    cluster_label <- meta$res.clusters@label
    scatter_weight <- meta$res.scatter@w
    cell_events <- meta$dat.clusters$clsEvents
    
    for( s in scatter ) {
# all cell-clusters in scatter-cluster s
        label <- sort(unique(cluster_label[scatter_label==s]))
        for( l in label ) {
# all scatter-labels for cell-clusters in meta-cluster l 
# t <- meta$res.scatter@label[meta$res.clusters@label==l]
            
            t <- unique(scatter_label[cluster_label==l])
            
            max_events <- 0
            max_scatter <- 0
            for( j in t ) {
                num_events <- scatter_weight[j] * 
                sum(cell_events[cluster_label==l & scatter_label==j])
                if( num_events > max_events ){
                    max_events <- num_events
                    max_scatter <- j
                }
            }
            if( max_scatter == s ) {
                ret <- c(ret,l)
            }
        }
    }
    ret
}


################################################################################

