## immunoClust accessors
setMethod("nobs", signature(object="immunoClust"),
function(object) {
    length(object@label)
})
setMethod("npar", signature(object="immunoClust"),
function(object) {
    length(object@parameters)
})
setMethod("ncls", signature(object="immunoClust"),
function(object, ...) {
    object@K
})
setMethod("weights", signature(object="immunoClust"),
function(object,cls=seq_len(ncls(object))) {
    object@w[cls]
})
setMethod("mu", signature(object="immunoClust"),
function(object, cls=seq_len(ncls(object)), par=seq_len(npar(object))) {
    object@mu[cls,par]
})
setMethod(sigma, signature(object="immunoClust"),
function(object, cls=seq_len(ncls(object)), par=seq_len(npar(object))) {
    object@sigma[cls,par,par]
})
setMethod("parameters", signature(object="immunoClust"),
function(object) {
    object@parameters
})
setReplaceMethod("parameters", 
signature=signature(object="immunoClust", value="character"),
function(object,value) {
    if( length(value) != ncol(object@mu))
    stop("length of value array does not match to number of parameter")
    object@parameters <- as.character(value)
    object
})
setMethod("label", signature(object="immunoClust"),
function(object) {
    object@label
})
setMethod("posterior", signature(object="immunoClust"),
function(object) {
    object@z
})
setMethod("events", signature(object="immunoClust"),
function(object,cls=seq_len(ncls(object)) ) {
    #ret <- sapply(cls, function(k) sum(!is.na(object@label) & object@label==k))
    ret <- vapply(cls, function(k) 
            sum(!is.na(object@label) & object@label==k), 0 )
    names(ret) <- sprintf("cls-%d",cls)
    ret
})

## immunoClust accessors
setMethod("subset", signature(x="immunoClust"),
function(x, cls=seq_len(ncls(x)), par=seq_len(npar(x)))
{
    P <- length(par)
    L <- length(cls)
    label <- rep(0, length(x@label))
    
    for( l in seq_len(L) ) {
        k <- cls[l]
        label[x@label==k] <- l
    }
    
    mu <- x@mu[cls,par]
    sigma <- x@sigma[cls,par,par]
    dim(mu) <- c(L,P)
    dim(sigma) <- c(L, P, P)
    
    y <- new("immunoClust", expName=x@expName,parameters=x@parameters[par],
        K=L,P=P,N=x@N,w=x@w[cls],mu=mu,sigma=sigma,
        z=matrix(0,nrow=0,ncol=0), label=label,
        logLike=x@logLike, BIC=x@BIC, ICL=x@ICL,
        history=x@history, state=x@state
        )
    
    desc <- attr(x, "desc")
    if( !is.null(desc) ) {
        attr(y,"desc") <- desc[par]
    }
    attr(y,"fcsName") <- attr(x,"fcsName")
    attr(y,"trans.a") <- attr(x,"trans.a")[par]
    attr(y,"trans.b") <- attr(x,"trans.b")[par]
    
    y
})

setMethod("transformParams", signature(object="immunoClust"),
function(object, scale=c(), offset=c())
{
    P <- npar(object)
    K <- ncls(object)
    if( length(scale) < P ) {
        scale <- c(scale, rep(1, P-length(scale)))
    }
    if( length(offset) < P ) {
        offset <- c(offset, rep(0, P-length(offset)))
    }
    
    y <- object
    for( k in seq_len(K))
    y@mu[k,] <- scale * object@mu[k,] + offset

    y
})
