.exprs <- function(x, varNames=NULL)
{
    if(length(varNames) == 0 ) {
        varNames <- colnames(x)
        varNames <- varNames[ varNames != "Time" ]
    }
    if (is(x, "flowFrame")) {
        if (length(varNames)==0) {
            y <- exprs(x)
            varNames <- colnames(y)
        }
        else {
            y <- as.matrix(exprs(x)[, varNames])
        }
    }
    else if (is(x, "matrix")) {
        if (length(varNames)==0) {
            y <- x
            if (length(colnames(x))==0) 
            varNames <- "Not Available"
            else varNames <- colnames(x)
        }
        else {
            y <- as.matrix(x[, varNames])
            if( ncol(y) == 1 ) {
                colnames(y) <- varNames
            }            
        }
    }
    else if (is(x, "data.frame")) {
        if (length(varNames)==0) {
            y <- as.matrix(x)
            varNames <- colnames(x)
        }
        else {
            y <- as.matrix(x[, varNames])
        }
    }
    else if (is(x, "vector")) {
        y <- matrix(x)
        if (length(varNames)==0) 
        varNames <- "Not Available"
    }
    else {
        stop(paste("Object ", as.character(x), 
                " is not of class flowFrame / matrix / data frame!"))
    }
    
    y
}

###
## meta.exp.data
##  collect data of all cell-clusters
###
meta.exprs <- function(exp, sub=c())
{
    N <- length(exp)
    
    params <- attributes(exp[[1]])$param
    if( is.null(params) ) {
        params <- exp[[1]]@varNames
    }
    
    if( is.null(sub) ) {
        sub <- seq_len(length(params))
        varNames <- params
    }
    else {
        varNames <- params[sub]
    }
    
    varDesc <- varNames
    
    map <- match(varNames, params)
    desc <- attr(exp[[1]], "desc")
    P <- length(varNames)
    K <- array(0, N)
    
    for( i in seq_len(N) ) {
        K[i] <- exp[[i]]@K
    }
    
    totK <- sum(K)
    
    W <- array(NA, dim=totK)
    M <- matrix(NA, nrow=totK, ncol=P)
    S <- matrix(NA, nrow=totK, ncol=P*P)
    cls <- rep("", totK)
    clsEvents <- rep(0, totK)
    expEvents <- rep(0, N)
    expNames <- rep("", N)
    removedEvents <- rep(0, N)
    if( !is.null(desc) ) {
        colnames(M) <- paste(sep="\n", varNames, desc[map])
        varDesc <- desc[map]
    }
    else {
        colnames(M) <- varNames
    }
    
    k <- 0
    for( i in seq_len(N) ) {
        res <- exp[[i]]
        expEvents[i] <- sum(!is.na(res@label))
        expNames[i] <- attr(res, "expName")
        removedEvents[i] <- sum(is.na(res@label))
        l <- k+1
        for( j in seq_len(K[i]) ) {
            clsEvents[k+j] <- sum(!is.na(res@label) & res@label==j)
        }
        
        k <- k + K[i]
        W[l:k] <- res@w
        M[l:k,] <- res@mu[,map]
        S[l:k,] <- res@sigma[,map,map]
        
        cls[l:k] <- paste(sep="", attr(exp[[i]], "expName"), "_", 1:K[i])
    }
    
    rownames(M) <- cls
    
    list("P"=P, "N"=N, "K"=K, "W"=W, "M"=M, "S"=S, "expNames"=expNames, 
        "expEvents"=expEvents, "removedEvents"=removedEvents, 
        "clsEvents"=clsEvents, "desc"=varDesc)
}


####
##  fixes compensation routine in flowCore
##  meanwhile obsolet
###
.cell.compensate <- function(x, parameters=NULL)
{
    spill <- spillover(x)$SPILL
    cols <- colnames(spill)
    if( !is.null(parameters) ) {
        sel <- cols %in% parameters
        spill[which(!sel),] <- 0
        spill[,which(!sel)] <- 0
        diag(spill) <- 1
    }
    e <- exprs(x)
    e[, cols] <- t(solve(t(spill))%*%t(e[,cols]))
    
    exprs(x) <- e
    x
}
### cell.compensate

removed.above <- function(fcs, parameters=NULL, N=NULL, max.count=10, max=NULL)
{
    
    dat <- fcs
## restrict number of events?
    if( !is.null(N) && N < nrow(dat) ) {
        dat <- dat[seq_len(N)]
    }
    else {
        N <- nrow(dat)
    }
    
    if( is.null(parameters) ) {
        parameters <- colnames(dat)
        parameters <- parameters[ parameters != "Time" ]
    }
    else {
        parameters <- parameters[!is.na(match(parameters,colnames(fcs)))]
    }
    
    y <- .exprs(dat, parameters) 
    
    removed <- matrix(0, ncol=ncol(y)+2, nrow=2)
    removed.not <- matrix(FALSE, ncol=ncol(y), nrow=nrow(y))
    
    rm.above <- rep(FALSE, nrow(y))
    if (max.count > -1) {
        if (is.null(max)[1]) 
        max <- apply(y, 2, max)
        for (p in seq_len(ncol(y)))
        if (sum(y[,p]>=max[p]) >= max.count) {
            removed[1,p] <-  sum(y[,p] >= max[p])
            rm.above <- rm.above | (y[,p] >= max[p])
            for( q in seq_len(ncol(y)) )
            if( q != p ) {
                removed.not[,q] <- removed.not[,q] | (y[,p] >= max[p])
            }
        }
        rm.sum <- sum(rm.above)
        for( p in seq_len(ncol(y)) ) {
            removed[2,p] <- sum(rm.above & !removed.not[,p])
        }
        
        removed[1,ncol(y)+1] <- rm.sum
        removed[1,ncol(y)+2] <- 100*rm.sum/nrow(y)
        
        removed[2,ncol(y)+1] <- sum(removed[2,seq_len(ncol(y))])
        removed[2,ncol(y)+2] <- 100*sum(removed[2,seq_len(ncol(y))])/nrow(y)
        
    }
        
    colnames(removed) <- c(parameters, "sum", "per ttl")
    rownames(removed) <- c("above", "above.only")
    
    removed
}
removed.below <- function(fcs, parameters=NULL, N=NULL, min.count=10, min=NULL)
{
    
    dat <- fcs
## restrict number of events?
    if( !is.null(N) && N < nrow(dat) ) {
        dat <- dat[seq_len(N)]
    }
    else {
        N <- nrow(dat)
    }
    
    if( is.null(parameters) ) {
        parameters <- colnames(dat)
        parameters <- parameters[ parameters != "Time" ]
    }
    else {
        parameters <- parameters[!is.na(match(parameters,colnames(fcs)))]
    }
    
    y <- .exprs(dat, parameters) 
    
    removed <- matrix(0, ncol=ncol(y)+2, nrow=2)
    removed.not <- matrix(FALSE, ncol=ncol(y), nrow=nrow(y))
    
    rm.below <- rep(FALSE, nrow(y))
    if (min.count > -1) {
        if (is.null(min)[1]) 
        min <- apply(y, 2, min)
        for (p in seq_len(ncol(y)))
        if (sum(y[,p]<=min[p]) >= min.count) {
            removed[1,p] <-  sum(y[,p] <= min[p])
            rm.below <- rm.below | (y[,p] <= min[p])
            for( q in seq_len(ncol(y)) ) 
            if( q != p ) {
                removed.not[,q] <- removed.not[,q] | (y[,p] <= min[p])
            }
        }
        rm.sum <- sum(rm.below)
        for( p in seq_len(ncol(y)) ) {
            removed[2,p] <- sum(rm.below & !removed.not[,p])
        }
        
        removed[1,ncol(y)+1] <- rm.sum
        removed[1,ncol(y)+2] <- 100*rm.sum/nrow(y)
        
        removed[2,ncol(y)+1] <- sum(removed[2,seq_len(ncol(y))])
        removed[2,ncol(y)+2] <- 100*sum(removed[2,seq_len(ncol(y))])/nrow(y)
        
    }
    
    
    colnames(removed) <- c(parameters, "sum", "per ttl")
    rownames(removed) <- c("below", "below.only")
    
    removed
}

