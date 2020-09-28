### 
#
#    cell (event) clustering
#
###

cell.process <- function(
fcs, parameters=NULL, apply.compensation=FALSE, classify.all=FALSE, 
N=NULL, min.count=10, max.count=10, min=NULL, max=NULL,
I.buildup=6, I.final=4, I.trans=I.buildup, 
modelName="mvt", tol=1e-5, bias=0.3,
sub.tol= 1e-4, sub.bias=bias, sub.thres=bias, sub.samples=1500, 
sub.extract=0.8, sub.weights=1, sub.standardize=TRUE,
trans.estimate=TRUE, trans.minclust=10, trans.a=0.01, trans.b=0.0, 
trans.parameters=NULL 
) {
## read file
##    fcs <- read.FCS(fcs_file)
    
    dat <- fcs
## restrict number of events?
    if( !is.null(N) && N < nrow(dat) ) {
        dat <- dat[seq_len(N)]
    }
    else {
        N <- nrow(dat)
    }
    
## filter over/under exposed
    if( is.null(parameters) ) {
        parameters <- colnames(dat)
        parameters <- parameters[ parameters != "Time" ]
    }
    else {
        parameters <- parameters[!is.na(match(parameters,colnames(fcs)))]
    }
    
    y <- .exprs(dat, parameters)
    rm.max <- rm.min <- rep(FALSE, nrow(y))
    if (max.count > -1) {
        if (is.null(max)[1]) 
        max <- apply(y, 2, max)
        for (p in seq_len(ncol(y)))  if (sum(y[,p]>=max[p]) >= max.count) 
        rm.max <- rm.max | (y[,p] >= max[p])
    }
    if (min.count > -1) {
        if (is.null(min)[1]) 
        min <- apply(y, 2, min)
        for (p in seq_len(ncol(y)))  if (sum(y[,p]<=min[p]) >= min.count) 
        rm.min <- rm.min | (y[,p] <= min[p])
    }
    inc <- !rm.max & !rm.min
    
    message("filtered from above:", sum(rm.max))
    message("filtered from below:", sum(rm.min), "\n")
    
    dat <- dat[inc]
    
    
## apply compensation    
    if( apply.compensation ) {
        dat <- .cell.compensate(dat)
    }
    
##     store fluorescence parmater information in description
    if( !is.null(trans.parameters) ) {
        trans.use <- !is.na(match(trans.parameters,colnames(fcs)))
        trans.parameters <- trans.parameters[trans.use]
        
        npar <- ncol(dat)
        id <- paste("P", seq_len(npar), "DISPLAY", sep="")
        #description(dat)[id] <- "LIN"
        keyword(dat)[id] <- "LIN"
        par <- match(trans.parameters, colnames(dat))
        id <- paste("P", par, "DISPLAY", sep="")
        #description(dat)[id] <- "LOG"
        keyword(dat)[id] <- "LOG"
    }
    
##
    if( trans.estimate ) {
        res <- cell.MajorIterationTrans(dat, parameters=parameters, 
                                    I.buildup=I.buildup, I.final=I.final, 
                                    I.trans=I.trans, modelName=modelName, 
                                    tol=tol, bias=bias, 
                                    sub.bias=sub.bias, sub.thres=sub.thres, 
                                    sub.tol=sub.tol, sub.samples=sub.samples, 
                                    sub.extract=sub.extract, 
                                    sub.weights=sub.weights, 
                                    sub.standardize=sub.standardize,
                                    trans.minclust=trans.minclust, 
                                    trans.a=trans.a)
    }
    else {
        m <- cell.InitialModel(dat, parameters=parameters,
                                trans.a=trans.a, trans.b=trans.b)
                                
        dat <- trans.ApplyToData(m, dat)
        
        res <- cell.MajorIterationLoop(dat, parameters=parameters, 
                                    I.buildup=I.buildup, I.final=I.final, 
                                    modelName=modelName, 
                                    tol=tol, bias=bias, 
                                    sub.bias=sub.bias, sub.thres=sub.thres, 
                                    sub.tol=sub.tol, sub.samples=sub.samples, 
                                    sub.extract=sub.extract, 
                                    sub.weights=sub.weights, 
                                    sub.standardize=sub.standardize)
        
## 
        attr(res, "trans.a") <- attr(m,"trans.a")
        attr(res, "trans.b") <- attr(m,"trans.b")
        attr(res, "trans.decade") <- attr(m, "trans.decade")
        attr(res, "trans.scale") <- attr(m, "trans.scale")
        
    }
## classify all events 
    if( classify.all ) {
        res@z <- matrix()
        dat <- fcs
        if( N < nrow(dat) ) {
            dat <- dat[seq_len(N)]
        }
        if( apply.compensation ) {
            dat <- .cell.compensate(dat)
        }
        t_dat <- trans.ApplyToData(res, dat)
        res <- cell.Classify(res, t_dat, modelName=modelName)
    }
    else {
        if( length(res@label) < N ) {
            label <- rep(NA,N)
            label[inc] <- res@label
            attr(res, "label") <- label
            attr(res, "removed.below") <- sum(rm.min)
            attr(res, "removed.above") <- sum(rm.max)
        }
    }
    
    attr(res, "bias") <- bias
    attr(res, "call") <- match.call()
    #attr(res, "fcsName") <- description(dat)$`FILENAME`
    attr(res, "fcsName") <- keyword(dat)$`FILENAME`
    
    par <- parameters(dat)
    inc <- match(res@parameters, par@data[,'name'])
    attr(res, "desc") <- par@data[inc, 'desc']
    
    res
}


### cell.InitialModel
cell.InitialModel <- function( 
dat, parameters=NULL, trans.a = 0.01, trans.b = 0.0, 
trans.decade=-1, trans.scale=1.0
) {
    
    if( is.null(parameters) ) {
        parameters <- colnames(dat)
    }
    
    N <- nrow(dat)
    P <- length(parameters)
    
    npar <- ncol(dat)
    
    desc <- attr(dat, "description")
    if( !is.null(desc) ) {
        par <- match(parameters, colnames(dat))
        
        id <- paste("P",par,"DISPLAY",sep="")
        display <- desc[id]
        
        sel <- display == "LOG"
        
        a <- rep(0.0,P)
        a[sel] <- trans.a
        b <- rep(0.0,P)
        b[sel] <- trans.b
    }
    else {
        a <- rep(0.0, P)
        a[] <- trans.a
        b <- rep(0.0, P)
        b[] <- trans.b
    }
    
    K <- 1
    z <- matrix(1, N, K)
    label <- rep(1, N)
    
# create model object 
    result <- new("immunoClust", expName="Initial Model", parameters=parameters,
                trans.a=a, trans.b=b, 
                trans.decade=trans.decade, trans.scale=trans.scale,
                K=1, N=N, P=P, z=z, label=label, 
                history=paste(1), state=rep(0,1)) 
    
    result
    
}
### cell.InitialModel

####


##
###
##    process (continue processing) clustering
##    data are transformed suitable before
###
cell.MajorIterationLoop <- function( 
dat, x=NULL, parameters=NULL, I.buildup=6, I.final=4, 
modelName="mvt", tol=1e-5, bias=0.3,
sub.bias=bias, sub.thres=0.0, sub.tol=1e-4, sub.samples=1500, 
sub.extract=0.8, sub.weights=1, sub.EM="MEt", 
sub.standardize=TRUE #, seed=1 
) {
    
    #set.seed(seed)
    s <- strptime(date(), "%a %b %d %H:%M:%S %Y")
    
    data <- dat
    N <- nrow(data)
    L <- I.buildup
    if( L > 0 ) {
        sam_t <- sample( seq_len(N), N/(2^L) )
        dat_t <- data[sam_t,]
        sam_t <- NULL
    }
    else {
        dat_t <- dat
    }
    
    if( is.null(x) ) {
## start with 1 cluster model
        res <- cell.ClustData(dat_t, 1, parameters=parameters, 
                            tol=tol, modelName=modelName)
    }
    else {
        res <- cell.Classify(x, dat_t)
        
    }
    
    #set.seed(seed)
    
    model <- NULL
    for( i in seq_len(I.buildup+I.final) ) {
        
        model <- cell.SubClustering(res, dat_t, thres=sub.thres, bias=sub.bias, 
                            B=100, tol=sub.tol, 
                            sample.weights=sub.weights, 
                            sample.EM=sub.EM,
                            sample.number=sub.samples, 
                            sample.standardize=sub.standardize, 
                            extract.thres=sub.extract, 
                            modelName=modelName)
        
        L <- L-1
        finished <- (L < (-1))
        if( L > 0 ) {
            sam_t <- sample( seq_len(N), N/(2^L) )
            dat_t <- data[sam_t,]
            sam_t <- NULL
        }
        else {
            dat_t <- data
        }
        
        if( model@K != res@K )
        finished <- FALSE
        
        message("Fit Model ", i, " of (", I.buildup, "/", I.buildup+I.final, 
                ") K=", res@K, "->", model@K, " N=", nrow(dat_t) )
        
        res <- NULL
        gc(verbose=FALSE, reset=TRUE)
        
        res <- cell.FitModel(model, dat_t, 
                            B=100, tol=tol, bias=bias, modelName=modelName)
        
        if( model@K != res@K )
        finished <- FALSE
        
        if( finished ) {
            message("Process completed at", i, "\n")
            break
        }
    }
    
    
    e <- strptime(date(), "%a %b %d %H:%M:%S %Y")
    message("Process ", I.buildup, "/", I.final, " takes ", 
            format(difftime(e,s,units="min"), digits=2), " minutes")
    res
}

cell.MajorIterationTrans <- function(
fcs, x=NULL, parameters=NULL, I.buildup=6, I.final=4, I.trans=I.buildup, 
modelName="mvt", tol=1e-5, bias=0.3,
sub.bias=bias, sub.thres=0.0, sub.tol=1e-4, sub.samples=1500, 
sub.extract=0.8, sub.weights=1, sub.EM="MEt", sub.standardize=TRUE, #seed=1,   
trans.minclust=5, trans.a=0.01, trans.decade=-1, trans.scale=1.0, 
trans.proc="vsHtransAw" 
) {
    
    data <- fcs
    s <- strptime(date(), "%a %b %d %H:%M:%S %Y")
    
    #set.seed(seed)
    N <- nrow(data)
    L <- I.buildup
    if( L > 0 ) {
        sam_t <- sample( seq_len(N), N/(2^L) )
        dat_t <- data[sam_t,]
        sam_t <- NULL
    }
    else {
        dat_t <- data
    }
    
    if( is.null(x) ) {
        
        model <- cell.InitialModel(dat_t, parameters=parameters, 
                                trans.a=trans.a, trans.decade=trans.decade, 
                                trans.scale=trans.scale )
        t_dat <- trans.ApplyToData(model, dat_t)
        res <- cell.ClustData(t_dat, 1, parameters=parameters, 
                            tol=tol, modelName=modelName)
        attr(res, "trans.a") <- attr(model, "trans.a")
        attr(res, "trans.b") <- attr(model, "trans.b")
        attr(res, "trans.decade") <- attr(model, "trans.decade")
        attr(res, "trans.scale") <- attr(model, "trans.scale")
        
    }
    else {
        t_dat <- trans.ApplyToData(x, dat_t)
        res <- cell.Classify(x, t_dat, modelName=modelName)
        
    }
    
    #set.seed(seed)
    
    model <- NULL
    for( i in seq_len(I.buildup+I.final) ) {
        L <- L-1
        finished <- (L < (-1))
        
        if( (i<=I.trans) && (res@K > trans.minclust) ) {
            
            t <- trans.FitToData(res, dat_t, 
                            B=10, certainty=0.3, proc=trans.proc)
            attr(res, "trans.a") <- t[1,]
            attr(res, "trans.b") <- t[2,]
            
        }
        
        model <- cell.SubClustering(res, t_dat, thres=sub.thres, bias=sub.bias, 
                                B=100, tol=sub.tol,
                                sample.weights=sub.weights,sample.EM=sub.EM, 
                                sample.number=sub.samples, 
                                sample.standardize=sub.standardize, 
                                extract.thres=sub.extract, 
                                modelName=modelName)
        
        if( L > 0 ) {
            sam_t <- sample( seq_len(N), N/(2^L) )
            dat_t <- data[sam_t,]
            sam_t <- NULL
        }
        else {
            dat_t <- data
        }
        
        t_dat <- trans.ApplyToData(res, dat_t)
        
        if( model@K != res@K )
        finished <- FALSE
        
        message("Fit Model ", i, " of (", I.buildup, "/", I.buildup+I.final, 
                ") K=", res@K, "->", model@K, " N=", nrow(t_dat) )
        res <- NULL
        gc(verbose=FALSE, reset=TRUE)
        
        res <- cell.FitModel(model, t_dat, B=50, tol=tol, bias=bias, 
                            modelName=modelName)
        
        if( model@K != res@K )
        finished <- FALSE
        if( finished ) {
            message("Process completed at", i, "\n")
            break
        }
    }
    
    e <- strptime(date(), "%a %b %d %H:%M:%S %Y")
    message("Major Iteration (Trans) (", I.buildup, "/", I.final, ") takes ", 
            format(difftime(e,s,units="min"), digits=2), " minutes")
    res
}


cell.classifyAll <- function(fcs, x, apply.compensation=FALSE)
{
##dat <- read.FCS(fcs_file)
    res <- x
    dat <- fcs
    N <- length(res@label)
    
    
    if( apply.compensation ) {
        dat <- .cell.compensate(dat)
    }
    if( N < nrow(dat) ) {
        dat <- dat[seq_len(N)]
    }
    t_dat <- trans.ApplyToData(res, dat)
    a_res <- cell.Classify(res, t_dat, modelName="mvt")
    
    #attr(a_res, "fcsName") <- description(dat)$`FILENAME`
    attr(a_res, "fcsName") <- keyword(dat)$`FILENAME`
    
    par <- parameters(dat)
    inc <- match(a_res@parameters, par@data['name'])
    attr(a_res, "desc") <- par@data[inc, 'desc']
    
    a_res
}




