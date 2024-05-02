.annotate.mergedClusters <- function(w, mu_, sigma_, cls)
{
    if( length(cls) == 0 ) {
        return( list("mu"=NULL, "sigma"=NULL) )
    }
    
    if( is.null(dim(mu_)) ) {
        P <- length(mu_)
        dim(mu_) <- c(1,P)
        dim(sigma_) <- c(1,P,P)
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

.annotate.buildModel <- function(pop, w, mu, sigma)
{
    r <- .annotate.mergedClusters(w, mu, sigma, pop$clusters)
    pop$M <- r$mu
    pop$S <- r$sigma
    if( !is.null(pop$childs) && length(pop$childs) > 0 ) {
        for( i in seq_len(length(pop$childs)) ) {
            pop$childs[[i]] <- .annotate.buildModel(pop$childs[[i]], 
                                    w, mu, sigma)
        }
    }
    pop
}


.annotate.newPop <- function(clusters=c(), childs=vector("list"),
    desc="new pop", partition=TRUE)
{
    list("clusters"=clusters, "childs"=childs, 
        "desc"=desc, "partition"=partition)
}

.annotate.setPop <- function(pop, pos, g, add=FALSE)
{
    if( is.null(pos) ) 
    return (g)
    
## get child classes not in here
    
    if( length(pop$clusters) == 0 ) {
        pop$clusters <- g$clusters
    }
    else {
        cls <- g$clusters[is.na(match(g$clusters,pop$clusters))]  
        if( length(cls) > 0 ) {
            message("found missing clusters in parent\n")
            pop$clusters <- sort(c(pop$clusters,cls))
        }
    }
    
    p <- pos[1]
    
    if( length(pop$childs) < p || is.null(pop$childs[[ p ]]) ) {
        if( !add ) 
        stop( "out of range in position level ", 
                paste(collapse=".", pop$pos), p, length(pop$childs), "\n")
        pop$childs[[ p ]] <- list("clusters"=c(), "childs"=vector("list"), 
                                "partition"=TRUE)
        pop$plot.childs <- TRUE
    }
    
    if( length(pos) == 1 ) {
        if( !is.null(g) ) {
            #message("add level ", p, " (", g$desc, ")\n\tto ", pop$desc, "\n")
            g$plot.subset <- pop$plot.subset
            g$parent.position <- pop$parent.position
            pop$childs[[ p ]] <- g
            pop$plot.childs <- TRUE
        }
        else {
            #message("remove level ", p, "\n\tfrom ", pop$desc, "\n")
            pop$childs[[ p ]] <- NULL
        }

    }
    else {
    
        pop$childs[[ p ]] <- .annotate.setPop(pop$childs[[p]], 
                                pos[2:length(pos)], g, add=add)
    }
    pop
}


.annotate.addPop <- function(pop, pos, clusters, desc)
{
    pop <- .annotate.setPop(pop, pos,.annotate.newPop(clusters, desc=desc), 
                            add=TRUE)
    pop <- .annotate.restructure(pop)
    pop
}


.annotate.getPop <- function(pop, pos)
{
    if( length(pos)==0  ) 
        return (pop)

    p <- pos[1]
    if( length(pos) > 1 ) {
        .annotate.getPop(pop$childs[[p]], pos[2:length(pos)])
    }
    else
    if( p > 0 && p <= length(pop$childs) ){
        pop$childs[[p]]
    }
    else {
        NULL
    }
    
}

.annotate.setParent <- function(pop, pos, parent=NULL, childs=TRUE)
{
    if( is.null(pop) ) {
        return(pop)
    }
    if( length(pos) == 0 ){ 
        if( is.null(parent) && childs ) {
            parent <- pop
        }
        else {
            if( is.numeric(parent) ) {
                pop$parent.desc <- NULL
                pop$parent.position <- parent
                pop$parent.clusters <- c()
            }
            else {
                pop$parent.desc <- parent$desc
                pop$parent.position <- parent$position
                pop$parent.clusters <- parent$clusters
            }
        }
        if( childs && !is.null(parent) ) {
            ## hierarchy of parent and unclassified clusters not brilliant
            if( is.numeric(parent) ) {
                pop$unclassified.parent.desc <- NULL
                pop$unclassified.parent.position <- parent
                pop$unclassified.parent.clusters <- c()
            }
            else {
                pop$unclassified.parent.desc <- parent$desc
                pop$unclassified.parent.position <- parent$position
                pop$unclassified.parent.clusters <- parent$clusters
            }
        }
        if( childs && length(pop$childs) > 0 ) {
            for( i in seq_len(length(pop$childs)) ) {
                pop$childs[[i]] <- 
                    .annotate.setParent(pop$childs[[i]], pos, parent,
                                        childs=childs)
            }
        }
        return (pop)
    }
    
    p <- pos[1]
    
    if( length(pos) == 1 ) {
        pop$childs[[p]] <- 
            .annotate.setParent(pop$childs[[p]], NULL, parent,
                                childs=childs)
        
    }
    else {
        pop$childs[[p]] <- 
            .annotate.setParent(pop$childs[[p]], pos[2:length(pos)], parent,
                                childs=childs)
    }
    pop
    
}


.annotate.getProp <- function(pop, pos, which)
{
    pop <- .annotate.getPop(pop, pos)
    pop[[which]]
}

.annotate.setProp <- function(pop, pos, which, value=NULL,
    for.level=TRUE, for.sublevels=FALSE)
{
    if( is.null(pop) ) {
        return(pop)
    }
    if( length(pos) == 0 ){
        if( for.level ) {
            pop[[which]] <- value
        }
        if( for.sublevels ) {
            for( i in seq_len(length(pop$childs)) ) {
                pop$childs[[i]] <-
                .annotate.setProp(pop$childs[[i]], c(), which, value,
                    for.level=TRUE, for.sublevels=TRUE)
            }
        }
        return (pop)
    }
    
    p <- pos[1]
    if( p < 1 || p > length(pop$childs) ) {
        stop( "out of range in position level ", paste(collapse=".", pop$pos), 
            " (", p, "/",length(pop$childs), ")\n")
        return (pop)
    }
    
    if( length(pos) == 1 ) {
        pop$childs[[p]] <-
        .annotate.setProp(pop$childs[[p]], NULL, which, value,
            for.level=for.level, for.sublevels=for.sublevels)
        
    }
    else {
        pop$childs[[p]] <- 
        .annotate.setProp(pop$childs[[p]], pos[2:length(pos)], which, value,
            for.level=for.level, for.sublevels=for.sublevels)
    }
    pop
}


.annotate.addClusters <- function(pop, clusters, pos)
{
    
## get child classes not in here
    
    if( length(pop$clusters) == 0 ) {
        pop$clusters <- clusters
    }
    else {
        cls <- clusters[is.na(match(clusters,pop$clusters))]  
        if( length(cls) > 0 ) {
            pop$clusters <- sort(c(pop$clusters,cls))
        }
    }

    if( length(pos) == 0 ) {
        return (pop)
    }
    
    p <- pos[1]
    
    if( length(pop$childs) < p || is.null(pop$childs[[ p ]]) ) {
        stop( "out of range in position level\n")
    }
    
    if( length(pos) == 1 ) {
        pop$childs[[ p ]] <- .annotate.addClusters(pop$childs[[p]],
                                clusters, NULL)
    }
    else {
        pop$childs[[ p ]] <- .annotate.addClusters(pop$childs[[p]],
                                clusters, pos[2:length(pos)])
    }
    pop
}

.annotate.removeClusters <- function(pop,  clusters, pos=c(), renumber=c())
{
    if( is.null(pop) ) 
    return(pop)
    
    if( length(pos) == 0 ) {
        pop$clusters <- pop$clusters[is.na(match(pop$clusters, clusters))]
        if( !is.null(renumber) ){
            pop$clusters <- match(pop$clusters,renumber)
        }
        if( length(pop$childs) > 0 ) {
            
            for( i in seq_len(length(pop$childs)) ) {
                pop$childs[[i]] <- .annotate.removeClusters(pop$childs[[i]],
                                        clusters, pos=pos, renumber=renumber)
            }
        }
        
        return (pop) 
    }
    
    p <- pos[1]
    
    if(length(pos) > 1 ) {
        pos <- pos[2:length(pos)]
    }
    else {
        pos <- NULL
    }
    
    pop$childs[[p]] <- .annotate.removeClusters( pop$childs[[p]], 
                            clusters, pos=pos)
    
    return(pop)
}

.annotate.clearClusters <- function(pop)
{
    if( is.null(pop) ) 
        return (pop)
    pop$clusters <- c()
    pop$unclassified.parent.clusters <- c()
    
    if( length(pop$childs) > 0 ) {
    for( i in seq_len(length(pop$childs)) ) {
        pop$childs[[i]] <- .annotate.clearClusters(pop$childs[[i]])
    }
    }

    return (pop) 

}

.annotate.moveClusters <- function(pop, clusters, pos)
{
    pop <- .annotate.removeClusters(pop, clusters)
    pop <- .annotate.addClusters(pop, clusters, pos)
    pop
}

.annotate.retrieveClusters <- function(pop, pos)
{
    if( is.null(pop) ) 
        return(c())

    if( length(pos) == 0 ) {
        return (sort(pop$clusters))
    }
    
    p <- pos[1]
        
    if(length(pos) > 1 ) {
        pos <- pos[2:length(pos)]
    }
    else {
        pos <- NULL
    }
    if( p > length(pop$childs))
        return(c())
        
    .annotate.retrieveClusters( pop$childs[[p]], pos)
}

.annotate.retrieveClassified <- function(pop, pos)
{
    if( is.null(pop) ) 
    return(c())
    
    if( length(pos) == 0 ) {
        cls <- c()
        if( length(pop$childs) > 0 ) {
            for( i in seq_len(length(pop$childs)) ) {
                if( !is.null(pop$childs[[i]]) ) {
                    cls <- c(cls, pop$childs[[i]]$clusters)
                }
            }
        }
        if( length(cls) > 1 )
            cls <- sort(cls)
        return (cls)
    }
    
    p <- pos[1]
    
    if(length(pos) > 1 ) {
        pos <- pos[2:length(pos)]
    }
    else {
        pos <- NULL
    }
    
    .annotate.retrieveClassified( pop$childs[[p]], pos)
}

.annotate.retrieveUnclassified <- function(pop, pos)
{
    if( is.null(pop) ) 
    return(c())
    
    if( length(pos) == 0 ) {
        
        clusters <- .annotate.retrieveClassified(pop,c())
        rest <- pop$clusters[ is.na( match(pop$clusters, clusters) ) ]
        
        return (sort(rest))
        
    }

    p <- pos[1]
    
    if(length(pos) > 1 ) {
        pos <- pos[2:length(pos)]
    }
    else {
        pos <- NULL
    }
    
    .annotate.retrieveUnclassified( pop$childs[[p]], pos)
}


.pop.restructure <- function(root, pop=root, pos=c(), subset=c(), 
    remove.empty=FALSE, depth=-1)
{
    if( is.null(pop) )
        return (pop)
    
    pop$position <- pos
    
    if( !is.null(subset) ) {
        plot.subset <- c()
        for( p in pop$plot.subset) {
            plot.subset <- c(plot.subset,which(subset==p))
        }
        if( length(plot.subset) <= 1 )
            plot.subset <- c()
        pop$plot.subset <- plot.subset
    }
    
    parent <- .annotate.getPop(root, pop$parent.position)
    pop$parent.desc <- parent$desc
    pop$parent.clusters <- parent$clusters
    pop$parent <- NULL
    
    if( !is.null(pop$unclassified.parent.position) ) {
    parent <- .annotate.getPop(root, pop$unclassified.parent.position)
    pop$unclassified.parent.desc <- parent$desc
    pop$unclassified.parent.clusters <- parent$clusters
    }
    
    if( length(pop$childs) > 0 && depth != 0 ) {
        for( i in seq_len(length(pop$childs)) ) {
            pop$childs[[i]] <- .pop.restructure(root, pop$childs[[i]], 
                                    c(pos,i), subset=subset,
            remove.empty=remove.empty, depth=depth-1)
        }
        if( remove.empty ) {
            for( i in length(pop$childs):1) {
                if( length(pop$childs[[i]]$clusters) == 0 )
                    pop$childs[[i]] <- NULL
            }
        }
    }
    else {
        pop$childs <- NULL
    }

    pop
}

.annotate.restructure <- function(root, remove.empty=FALSE, depth=-1)
{
    root <- .pop.restructure(root, remove.empty=remove.empty, depth=depth)
    root$parent.desc <- NULL
    root$parent.clusters <- NULL
    
    root
}


.annotate.clustering <-
function(meta.src, meta.dst, param.subset=c(),
meta.iter=10, meta.bias=0.1,
meta.alpha=0.5, norm.method=3, norm.blur=2, norm.minG=10, verbose=FALSE)
{
    anno.src <- meta.src$gating
    anno.dst <- .annotate.clearClusters(anno.src)
    anno.dst$clusters <- seq_len(meta.dst$res.clusters@K)
    
    ## co-clustering and normalization
    anno.exp <- vector("list", 2)
    anno.exp[[1]] <- meta.src$res.clusters
    anno.exp[[2]] <- meta.dst$res.clusters
    
    dat <- meta.exprs(anno.exp)
    
    ## clsEvents
    l <- 0
    exp.res <- meta.src$res.clusters
    exp.dat <- meta.src$dat.clusters
    for( k in seq_len(exp.res@K) ) {
        cls <- which( exp.res@label == k )
        dat$clsEvents[l+k] <- sum(exp.dat$clsEvents[cls])
    }
    l <- l + exp.res@K
    exp.res <- meta.dst$res.clusters
    exp.dat <- meta.dst$dat.clusters
    for( k in seq_len(exp.res@K) ) {
        cls <- which( exp.res@label == k )
        dat$clsEvents[l+k] <- sum(exp.dat$clsEvents[cls])
    }
    
    
    ##
    ## doclustering
    res <- meta.Clustering(dat$P, dat$N, dat$K,
    dat$clsEvents, dat$M, dat$S,
    bias=meta.bias, I.iter=meta.iter, EM.method=20,
    alpha=meta.alpha, norm.method=norm.method,
    norm.blur=norm.blur, norm.minG=norm.minG)
    ## do normalization
    dat.norm <- meta.Normalize(dat$P, dat$N, dat$K, dat$clsEvents,
    dat$M, dat$S, res@K, res@z,
    method=norm.method)
    
    
    ## build normalized model
    src.M <- dat.norm$M[seq_len(dat$K[1]),]
    src.S <- dat.norm$S[seq_len(dat$K[1]),]
    src.W <- dat.norm$W[seq_len(dat$K[1])]
    
    dst.M <- dat.norm$M[seq_len(dat$K[2])+dat$K[1],]
    dst.S <- dat.norm$S[seq_len(dat$K[2])+dat$K[1],]
    dst.W <- dat.norm$W[seq_len(dat$K[2])+dat$K[1]]
    
    P <- ncol(src.M)
    dim(src.S) <- c(dat$K[1],P,P)
    dim(dst.S) <- c(dat$K[2],P,P)
    
    anno.src <- .annotate.buildModel(anno.src, src.W, src.M, src.S)
    
    find.pop <- function(anno, M, S, tst) {
        if( length(anno$clusters)==0 )
        return(tst)
        
        ret <- tst
        t <- tst
        t$prop <- bhattacharyya.prob(anno$M, anno$S, M, S)
        if( !is.na(t$prop) && t$prop > tst$prop ) {
            ret <- t
        }
        else{
        }
        if( length(anno$childs) ) {
            clusters <- c()
            for( i in seq_len(length(anno$childs)) ) {
                t <- find.pop(anno$childs[[i]], M, S,
                list("pos"=c(tst$pos,i), "prop"=ret$prop))
                if( t$prop > ret$prop ) {
                    ret <- t
                }
                clusters <- c(clusters, anno$childs[[i]]$clusters)
            }
            ## not jet tested
            rest <- anno$clusters[ is.na( match(anno$clusters, clusters) ) ]
            
        }
        else {
            ## test all clusters in pop
            rest <- anno$clusters
        }
        ## test not jet respected meta-clusters
        t <- tst
        for( c in rest ) {
            t$prop <- bhattacharyya.prob(src.M[c,], src.S[c,,], M, S)
            if(t$prop > ret$prop ) {
                ret <- t
            }
        }
        
        ret
    }
    
    prop <- vector("list", dat$K[2])
    for( k in seq_len(dat$K[2]) ) {
        t <- find.pop(anno.src, dst.M[k,], dst.S[k,,], 
                        list("pos"=c(), "prop"=0))
        if( verbose )
        message("cls ", k, " moveto ", t$pos, t$prop)
        prop[[k]] <- t
        anno.dst <- .annotate.moveClusters(anno.dst, c(k), t$pos)
    }
    
    if( !is.null(param.subset) ) {
        src.M <- src.M[,param.subset]
        src.S <- src.S[,param.subset,param.subset]
        dst.M <- dst.M[,param.subset]
        dst.S <- dst.S[,param.subset,param.subset]
        
        anno.src <- .annotate.buildModel(anno.src, src.W, src.M, src.S)
        for( k in seq_len(dat$K[2]) ) {
            t <- find.pop(anno.src, dst.M[k,], dst.S[k,,], 
                            list("pos"=c(), "prop"=0))
            if( t$prop > prop[[k]]$prop ) {
                if( verbose )
                message("cls ", k, " move from ", prop[[k]]$pos, prop[[k]]$prop,
                " to ", t$pos, t$prop)
                if( (t$prop / prop[[k]]$prop) > 2 )
                anno.dst <- .annotate.moveClusters(anno.dst, c(k), t$pos)
            }
        }
        
    }
    
    anno.dst <- .annotate.restructure(anno.dst)
    
    meta.dst$gating <- anno.dst
    
    class(meta.dst) <- "immunoMeta"
    class(meta.dst$gating) <- "immunoAnno"
    
    meta.dst
    
}
