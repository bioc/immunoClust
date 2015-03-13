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
        sub <- 1:length(params)
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
    
    for( i in 1:N ) {
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
    if( !is.null(desc) ) {
        colnames(M) <- paste(sep="\n", varNames, desc[map])
        varDesc <- desc[map]
    }
    else {
        colnames(M) <- varNames
    }
    
    k <- 0
    for( i in 1:N ) {
        res <- exp[[i]]
        expEvents[i] <- sum(!is.na(res@label))
        expNames[i] <- attr(res, "expName")
        l <- k+1
        for( j in 1:(K[i]) ) {
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
        "expEvents"=expEvents, "clsEvents"=clsEvents, "desc"=varDesc)
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
        spill[which(!sel),which(!sel)] <- 1
    }
    e <- exprs(x)
    e[, cols] <- t(solve(t(spill))%*%t(e[,cols]))
    
    exprs(x) <- e
    x
}
### cell.compensate


