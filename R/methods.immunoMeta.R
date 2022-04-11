## immunoMeta accessors
setMethod("nsam", signature(object="immunoMeta"),
function(object) {
    object$dat.clusters$N
})
setMethod("sam_ncls", signature(object="immunoMeta"),
function(object, for.samples=seq_len(nsam(object))) {
    object$dat.clusters$K[for.samples]
})
setMethod("sam_clsWeights", signature(object="immunoMeta"),
function(object) {
    object$dat.clusters$W
})
setMethod("sam_clsMu", signature(object="immunoMeta"),
function(object) {
    object$dat.clusters$M
})
setMethod("sam_clsSigma", signature(object="immunoMeta"),
function(object) {
    object$dat.clusters$S
})
setMethod("sam_clsEvents", signature(object="immunoMeta"),
function(object) {
    object$dat.clusters$clsEvents
})
setMethod("nobs", signature(object="immunoMeta"),
function(object) {
    length(object$res.clusters@label)
})
setMethod("npar", signature(object="immunoMeta"),
function(object) {
    length(object$res.clusters@parameters)
})
setMethod("ncls", signature(object="immunoMeta"),
function(object) {
    object$res.clusters@K
})
setMethod("weights", signature(object="immunoMeta"),
function(object, cls=seq_len(ncls(object))) {
    ret <- object$res.clusters@w[cls]
    names(ret) <- sprintf("cls-%d",cls)
    ret
})
setMethod("mu", signature(object="immunoMeta"),
function(object, cls=seq_len(ncls(object)), par=seq_len(npar(object))) {
    ret <- object$res.clusters@mu[cls,par]
    #if( is.null(dim(ret)) ) {
    #    dim(ret) <- c(length(cls), length(par))
    #}
    if( !is.null(dim(ret) ) ) {
        rownames(ret) <- sprintf("cls-%d",cls)
        colnames(ret) <- sprintf("par-%d",par)
    }
    else {
        names(ret) <- sprintf("cls-%d",cls)
    }
    ret
})
setMethod("sigma", signature(object="immunoMeta"),
function(object, cls=seq_len(ncls(object)), par=seq_len(npar(object))) {
    ret <- object$res.clusters@sigma[cls,par,par]
    ret
})
setMethod("parameters", signature(object="immunoMeta"),
function(object) {
    object$res.clusters@parameters
})
setReplaceMethod("parameters", 
signature=signature(object="immunoMeta", value="character"),
function(object, value) {
    if( length(value) != ncol(object$res.clusters@mu))
    stop("length of value array does not match to number of parameter")
    
    object$res.clusters@parameters <- as.character(value)
    ## needs to fit
    colnames(object$dat.clusters$M) <- object$res.clusters@parameters
    object
})
setMethod("label", signature(object="immunoMeta"),
function(object, for.sample=NA) {
    if( is.na(for.sample) )
    return(object$res.clusters@label)
    
    #if( !is.integer(for.sample) )
    if( abs(for.sample-round(for.sample)) > .Machine$double.eps )
    stop("for.sample option has to be an integer")
    
    K <- object$dat.clusters$K
    
    if( for.sample < 1 || for.sample > length(K) )
    stop("for.sample option is out of range")
    
    sl <- sum(K[seq_len(for.sample-1)])+1
    el <- sum(K[seq_len(for.sample)])
    
    object$res.clusters@label[sl:el]
    
})
setMethod("aposteriori", signature(object="immunoMeta"),
function(object){
    object$res.clusters@z
})
setMethod("events", signature(object="immunoMeta"),
function(object,cls=seq_len(ncls(object)), for.sample=NA) {
    K <- object$dat.clusters$K
    N <- length(K)
    
    if( !is.na(for.sample) ) {
        if( for.sample < 1 || for.sample > length(K) )
        stop("for.sample option is out of range")
        
        n <- for.sample
        
        ret <- c()
        for( j in cls ) {
            k <- which(label(object)==j)
            inc <- sum(K[seq_len(n-1)]) < k & k <= sum(K[seq_len(n)])
            ret <- c(ret, sum(object$dat.clusters$clsEvents[ k[inc] ]))
        }
        names(ret) <- sprintf("cls-%d",cls)
        
        return (ret)
    }
    
    ret <- matrix(0,nrow=N,ncol=0)
    
    for( j in cls ) {
        k <- which(label(object)==j)
        col <- c()
        for( n in seq_len(N) ) {
            inc <- sum(K[seq_len(n-1)]) < k & k <= sum(K[seq_len(n)])
            col <- c(col,sum(object$dat.clusters$clsEvents[ k[inc] ] ))
        }
        ret <- cbind(ret,col)
    }
    
    colnames(ret) <- sprintf("cls-%d",cls)
    rownames(ret) <- sprintf("exp-%d",seq_len(N))
    ret
})

setMethod("prop", signature(object="immunoMeta"),
function(object, name="", pos=c())
{
    pop <- .annotate.getPop(object$gating, pos)
    pop[[name]]
})

#setReplaceMethod("prop", signature=signature(object="immunoMeta", value="ANY"),
#definition=function(object, name="", pos=NA, for.level=TRUE,
#for.sublevels=FALSE, ..., value )
#{
#     if( !all(is.integer(pos)) ) {
#           stop("level position has to be specified as integer array")
#    }
#    object$gating <- .annotate.setProp(object$gating, pos, name, value,
#        for.level=for.level, for.sublevels=for.sublevels)
#    object
#})

"prop<-.immunoMeta" <-
function(object, name, pos, for.level=TRUE, for.sublevels=FALSE, ..., value)
{
    object$gating <- .annotate.setProp(object$gating, pos, name, value,
    for.level=for.level, for.sublevels=for.sublevels)
    object
}

setMethod("desc", signature(object="immunoMeta"),
function(object, pos=c())
{
    #   cat("immunoMeta", which, "\n")
    pop <- .annotate.getPop(object$gating, pos)
    pop[["desc"]]
})


"desc<-.immunoMeta" <-
function(object, pos, ..., value)
{
    object$gating <- .annotate.setProp(object$gating, pos, "desc", value)
    object
}

setMethod("descFull", signature(object="immunoMeta"),
function(object, pos=c())
{
    desc <- desc(object,c())
    len <- length(pos)
    for( i in seq_len(len) ) {
        desc <- c(desc,desc(object,pos[seq_len(i)]))
    }
    paste(collapse="_", desc)
})

setMethod("level", signature(object="immunoMeta"),
function(object, pos=c())
{
    pop <- .annotate.getPop(object$gating, pos)
    pop
})
#setReplaceMethod("level",signature=signature(object="immunoMeta",value="ANY"),
#definition=function(object, pos=NA, ..., value)
#{
#    if( !all(is.integer(pos)) ) {
#        stop("level position has to be specified as integer array")
#    }
#    object$gating <- .annotate.setPop(object$gating, pos, value)
#    object$gating <- .annotate.restructure(object$gating)
#    object
#})

"level<-.immunoMeta" <- function(object, pos, ..., value)
{
    object$gating <- .annotate.setPop(object$gating, pos, value)
    object$gating <- .annotate.restructure(object$gating)
    object
}

setMethod("findLevel", signature(object="immunoMeta"),
function(object, cls)
{
    .find.level <- function(pop, pos) {
        
        if( is.null(pop) )
        return(pos)
        
        if( !(cls %in% pop$clusters) )
        return(pos)
        
        #if( verbose )
        #message("search in ", paste(collapse=".", pop$position))
        
        pos <- pop$position
        
        if( length(pop$childs) > 0 ) {
            for( i in seq_len(length(pop$childs)) )
            pos <- .find.level(pop$childs[[i]], pos)
        }
        
        pos
    }
    
    .find.level(object$gating, c())
})

setMethod("clusters", signature(object="immunoMeta"),
function(object, pos)
{
    cls <- .annotate.retrieveClusters(object$gating, pos)
    cls
})
setMethod("classified", signature(object="immunoMeta"),
function(object, pos)
{
    cls <- .annotate.retrieveClassified(object$gating, pos)
    cls
})
setMethod("unclassified", signature(object="immunoMeta"),
function(object, pos)
{
    cls <- .annotate.retrieveUnclassified(object$gating, pos)
    cls
})
## immunoMeta accessors

## immunoMeta manipulators
#setReplaceMethod("addLevel",
#signature=signature(object="immunoMeta", value="ANY"),
#definition=function(object, pos=NA, desc="new level", ..., value)
#{
#    cat("addLvel", pos,"\n")
#    if( length(pos) > 0 && !all(is.integer(pos)) ) {
#        stop("level position has to be specified integer array")
#    }
#    if( !all(value %in% seq_len(ncls(object))) ) {
#        stop("some level clusters are not in meta cluster range")
#    }#
#
#    object$gating <- .annotate.addPop(object$gating, pos, value, desc)
#    object
#})

"addLevel<-.immunoMeta" <- function(object, pos, desc="new level", ..., value)
{
    if( !all(value %in% seq_len(ncls(object))) ) {
        stop("some level clusters are not in meta cluster range")
    }

    object$gating <- .annotate.addPop(object$gating, pos, value, desc)
    object
}

#setReplaceMethod("move",
#signature=signature(object="immunoMeta", value="integer"),
#definition=function(object, pos=NA, add=FALSE, ..., value)
#{
#    if( length(pos) > 0 && !all(is.integer(pos)) ) {
#        stop("level position has to be specified as integer array")
#    }
#
#    if( !all(value %in% seq_len(ncls(object))) ) {
#        stop("some clusters are not in cluster range")
#    }
#
#    if( add )
#    object$gating <- .annotate.addClusters(object$gating, value, pos)
#    else
#    object$gating <- .annotate.moveClusters(object$gating, value, pos)
#
#    object
#})

"move<-.immunoMeta" <- function(object, pos, add=FALSE, ..., value)
{
    if( !all(value %in% seq_len(ncls(object))) ) {
        stop("some clusters are not in cluster range")
    }

    if( add )
    object$gating <- .annotate.addClusters(object$gating, value, pos)
    else
    object$gating <- .annotate.moveClusters(object$gating, value, pos)

    object
}
#setReplaceMethod("remove",
#signature=signature(object="immunoMeta", value="ANY"),
#definition=function(object, pos=NA, ..., value)
#{
    #    if( !all(value %in% 1:x$res.clusters@K) ) {
    #       stop("some clusters are not in cluster range")
    #   }
#     if( length(pos) > 0 && !all(is.integer(pos)) ) {
#          stop("level position has to be specified as integer array")
#   }
#   if( value == "all" ) {
#       object$gating <- .annotate.clearClusters(object$gating)
#       object$gating <- .annotate.addClusters(object$gating,
#                           seq_len(ncls(object)), c())
#   }
#   else {
#       object$gating <- .annotate.removeClusters(object$gating, value, pos=pos)
#   }
#   object
#})

"remove<-.immunoMeta" <- function(object, pos, ..., value)
{
    #    if( !all(value %in% 1:x$res.clusters@K) ) {
    #       stop("some clusters are not in cluster range")
    #   }
    
    if( value == "all" ) {
        object$gating <- .annotate.clearClusters(object$gating)
        object$gating <- .annotate.addClusters(object$gating,
                            seq_len(ncls(object)), c())
    }
    else {
        object$gating <- .annotate.removeClusters(object$gating, value, pos=pos)
    }
    object
}

#setReplaceMethod("parent",
#signature=signature(object="immunoMeta", value="ANY"),
#definition=function(object, pos=NA, sub.levels=TRUE, ..., value)
#{
#    if( length(pos) > 0 && !all(is.integer(pos)) ) {
#           stop("level position has to be specified as integer array")
#    }
#    if( !is.numeric(value) &&
#          !all(value$clusters %in% seq_len(ncls(object)) ) ) {
#          stop("some parent clusters are not in cluster range")
#      }
#
#      object$gating <- .annotate.setParent(object$gating, pos, parent=value,
#                         childs=sub.levels)
#      object
#
#})
"parent<-.immunoMeta" <- function(object, pos, sub.levels=TRUE, ..., value)
{
    if( !is.numeric(value) &&
        !all(value$clusters %in% seq_len(ncls(object)) ) ) {
        stop("some parent clusters are not in cluster range")
    }

    object$gating <- .annotate.setParent(object$gating, pos, parent=value,
                        childs=sub.levels)
    object
}

setMethod("subset", signature(x="immunoMeta"),
function(x, cls=seq_len(ncls(x)), par=seq_len(npar(x)))
{
    y <- x
    if( length(cls) != ncls(x) ) {
        K <- ncls(x)
        cls <- sort(cls)
        rem <- rep(TRUE,K)
        rem[cls] <- FALSE
        
        P <- x$res.clusters@P
        w <- x$res.clusters@w[cls]
        clsEvents <- x$dat.clusters$clsEvents[cls]
        m <- x$res.clusters@mu[cls,]
        colnames(m) <- colnames(x$dat.clusters$M)
        s <- x$res.clusters@sigma
        dim(s) <- c(K,P*P)
        #cat("dim S", dim(s), "\n")
        s <- s[cls,]
        
        K <- length(cls)
        #cat("dim S", dim(s), "K", K, "P", P, "\n")
        y$dat.clusters$K <- K
        y$dat.clusters$W <- w
        y$dat.clusters$M <- m
        y$dat.clusters$S <- s
        y$dat.clusters$clsEvents <- clsEvents
        
        dim(s) <- c(K,P,P)
        y$res.clusters@K <- K
        y$res.clusters@w <- w
        y$res.clusters@mu <- m
        y$res.clusters@sigma <- s
        
        
        if( ncls(x) == nobs(x) ) {
            renumber <- which(!rem)
            y$res.clusters@label <- match(x$res.clusters@label[cls], renumber)
        
            y$gating <- .annotate.removeClusters(x$gating,which(rem), pos=c(),
                renumber=renumber)
        }
        else {
            y$gating <- .annotate.removeClusters(x$gating,which(rem), pos=c())
        }
        
    }
    if(length(par) != npar(x)) {
    params <- par
    ## restrict to subset parameter
    res <- y$res.clusters
    res@mu <- res@mu[,params]
    res@sigma <- res@sigma[,params,params]
    res@P <- length(params)
    #res@desc <- res@desc[params]
    res@parameters <- res@parameters[params]
    dat <- y$dat.clusters
    dat$M <- dat$M[, params]
    dim(dat$S) <- c(nrow(dat$S), dat$P, dat$P)
    dat$S <- dat$S[,params,params]
    dat$P <- length(params)
    dim(dat$S) <- c(nrow(dat$S),dat$P*dat$P)
    dat$desc <- dat$desc[params]
    
    sigma <- attr(res, "sigma.scaled")
    if( !is.null(sigma)){
        sigma <- sigma[,params,params]
        dim(sigma) <- c(res@K,res@P,res@P)
        attr(res,"sigma.scaled") <- sigma
    }
    
    y <- immunoMeta(res, dat, y$gating)
    y$gating <- .pop.restructure(y$gating, subset=params)
    #y$gating <- .annotate.buildModel(y$gating, y$res.clusters)
    y$gating <- .annotate.buildModel(y$gating, weights(y), mu(y), sigma(y))
    if( !is.null(prop(x,"pscales",c())) )
        prop(y,"pscales",c()) <- prop(x,"pscales",c())[params]
    }
    
    y
})

setMethod("transformParams", signature(object="immunoMeta"),
function(object, scale=c(), offset=c(), scale.sigma=FALSE)
{
    
    #   res <- x$res.clusters
    #   dat <- x$dat.clusters
    P <- object$dat.clusters$P
    K <- object$res.clusters@K
    totK <- sum(object$dat.clusters$K)
    if( length(scale) < P ) {
        scale <- c(scale, rep(1, P-length(scale)))
    }
    if( length(offset) < P ) {
        offset <- c(offset, rep(0, P-length(offset)))
    }
    
    y <- object
    for( k in seq_len(K) )
    y$res.clusters@mu[k,] <- scale * object$res.clusters@mu[k,] + offset
    for( k in seq_len(totK) )
    y$dat.clusters$M[k,] <- scale * object$dat.clusters$M[k,] + offset
    
    if( scale.sigma ) {
        s <- ((scale) %*% t((scale))) ## s[p,q] = scale[p] * scale[q]
        for( k in seq_len(K) )
        y$res.clusters@sigma[k,,] <- s * y$res.clusters@sigma[k,,]
    }
    
    y <- finalize(y)
    y
})

setMethod("finalize", signature(object="immunoMeta"),
function(object, remove.empty=FALSE, depth=-1)
{
    object$gating <- .annotate.restructure(object$gating,
                    remove.empty=remove.empty, depth=depth)
    object$gating <- .annotate.buildModel(object$gating,
        weights(object), mu(object), sigma(object))
    object
})

"transfer<-.immunoMeta" <- function(object, value)
{
    object <- .annotate.clustering(value, object)
    object <- finalize(object)
    object
}

#setReplaceMethod("transfer",
#signature=signature(object="immunoMeta",  value="immunoMeta"),
#definition=function(object, ..., value)
#{
    #object <- .annotate.clustering(value, object)
    #object <- finalize(object)
    #object
#    transfer(object) <- value
#})

setMethod("clusterDist", signature="immunoMeta",
function(object, cls, lvl, par=seq_len(npar(object)) ) {
    
    d <- vapply(cls, function(k) bhattacharyya.dist(prop(object,"M", lvl)[par],
        prop(object, "S", lvl)[par,par],
        mu(object,k,par), sigma(object,k,par)), 0 )
    
    names(d) <- sprintf("cls-%d",cls)
    d
})

setMethod("clusterCoeff", signature="immunoMeta",
function(object, cls, lvl, par=seq_len(npar(object)) ) {

    d <- vapply(cls, function(k) bhattacharyya.coeff(prop(object,"M", lvl)[par],
        prop(object, "S", lvl)[par,par],
        mu(object,k,par), sigma(object,k,par)), 0 )
    
    
    names(d) <- sprintf("cls-%d",cls)
    d
})

setMethod("clusterProb", signature="immunoMeta",
function(object, cls, lvl, par=seq_len(npar(object)) ) {
        
    d <- vapply(cls, function(k) bhattacharyya.prob(prop(object,"M", lvl)[par],
        prop(object, "S", lvl)[par,par],
        mu(object,k,par), sigma(object,k,par)), 0 )
    
    
    names(d) <- sprintf("cls-%d",cls)
    d
})

