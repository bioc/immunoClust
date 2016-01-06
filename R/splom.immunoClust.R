
.ellipse.panel.xyplot <- function(x, y, frame, subset, include, 
col=include+1, ecol=col, elty=1, eqtl=0.9, npoints=501, add=FALSE, ...)
{
    col <- matrix(col, length(include))
    
    panel.points( c(0,5,0,5), c(0,0,5,5), type="p", pch="+", col="gray" )
    
    
    elty <- matrix(elty, length(include))
    
    cc <- qt(eqtl, 5)
    
    j <- 0
    for (i in include) {
        j <- j+1
        eigenPair <- eigen(frame$sigma[i,subset,subset])
        l1 <- sqrt(eigenPair$values[1]) * cc
        l2 <- sqrt(eigenPair$values[2]) * cc
        angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1]) * 180/pi
        
        panel.points(.ellipsePoints(a=l1, b=l2, alpha=angle, 
                                    loc=frame$mu[i,subset], n=npoints), 
                    type="l", lty=elty[j], col=col[j])
        
    }
    
}

.ellipse.panel.splom <- function(x, y, frame, ...)
{
    cv <- current.viewport()$name
    dims <- as.numeric(strsplit(cv, ".", fixed=TRUE)[[1]][2:3])    
    
    .ellipse.panel.xyplot(x=x, y=y, frame=frame, subset=c(dims[1],dims[2]), ...)
}

.clust.prepanel <- function(x, y, frame, subset,...)
{
    range <- frame$range[subset]
    
    xlim <- c(range[1,1], range[2,1])
    ylim <- c(range[1,2], range[2,2])
    
    if( subset[1] > 2 ) {
        xlim <- c(-0.5, 5)
    }
    if( subset[2] > 2 ) {
        ylim <- c(-0.5, 5)
    }
    
    list("xlim"=xlim, "ylim"=ylim)
}


.clust.panel.xyplot <- function(x, y, frame, subset, include, 
col=include+1, pch=".",cex=0.6, 
ellipse=FALSE, ellipse.force=FALSE, ellipse.quantile=0.9, ecol=col, elty=1,
show.rm=FALSE, col.rm=1, pch.rm=1, cex.rm=0.6, 
npoints=501, add=FALSE, gates=NULL, lda=NULL, ...)
{
    
    label <- frame$label
    
    if( !is.null(frame$range) ) {
        range <- frame$range[subset]
        flagFiltered <- is.na(label) | (x < range[1,1]) | (y < range[1,2])
        panel.points(c(-0.5, 0,range[2,1],0,range[2,1]), 
                    c(-0.5, 0,0,range[2,2],range[2,2]), 
                    type="p", pch="+", col="gray")
    }
    else {
        flagFiltered <- is.na(label)
    }
    
    col <- matrix(col, length(include))
    pch <- matrix(pch, length(include))
    cex <- matrix(cex, length(include))
    
    
# if( !add ) plot(x, y, type="n", main=NULL, ...)
#    panel.points(x,y, type="n");
    
    if( !is.null(gates) && !ellipse.force ) {
        thres <- gates[subset,]
        for( j in 1:length(thres[1,]) ){
            if( !is.na(thres[1,j]) || !is.na(thres[2,]) )
            ellipse <- FALSE
        }
    }
    
    
# plot ellipses
    if (ellipse) {
        ecol <- matrix(ecol, length(include))
        elty <- matrix(elty, length(include))
        
        cc <- qt(ellipse.quantile, 5)
        
        j <- 0
        cc <- rep(cc, length.out=frame$K)
        for (i in include) {
            j <- j+1
            eigenPair <- eigen(frame$sigma[i,subset,subset])
            l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
            l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
            angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1]) 
            angle <- angle * 180/pi
            
            panel.points(.ellipsePoints(a=l1[i], b=l2[i], alpha=angle, 
                                        loc=frame$mu[i,subset], n=npoints), 
                        type="l", lty=elty[j], col=ecol[j])
            
        }  
    }
    
# plot gates    
    if( !is.null(gates) ) {
        x.limits = c(min(x[!flagFiltered],-1), max(x[!flagFiltered],10))
        y.limits = c(min(y[!flagFiltered],-1), max(y[!flagFiltered],10))
        thres <- gates[subset,]
        for( j in 1:length(thres[1,]) ) {   
            if( !is.na(thres[1,j]) ) {
                panel.points(c(thres[1,j],thres[1,j]), y.limits, 
                            type="l", col="red")
            }
        }
        for( j in 1:length(thres[2,]) ) {   
            if( !is.na(thres[2,j]) ) {
                panel.points(x.limits, c(thres[2,j],thres[2,j]), 
                            type="l", col="red")
            }
        }
    }
    
# plot LDA
    if( !is.null(lda) ) {
        l1 <- lda[[ subset[1] ]]
## l2 <- lda[[ subset[2] ]]
        if( !is.null(l1) ) {
            x.limits = c(min(x[!flagFiltered],-1), max(x[!flagFiltered],10))
            y.limits = c(min(y[!flagFiltered],-1), max(y[!flagFiltered],10))
            
            center <- l1$center[subset]
            vector <- l1$vector[subset]
            if( sum(is.na(vector))== 0 ) {
                x <- x.limits
                y <- c( center[2] + (x.limits[1]-center[1]) 
                        / vector[1] * vector[2],
                        center[2] + (x.limits[2]-center[1]) 
                        / vector[1] *vector[2] )
                panel.points(x.limits, y, type="l", col="green")
                
# center <- l1$center[subset]
# vector <- l1$vector[subset]
# y <- y.limits
# x <- c(center[1] - (y.limits[1]-center[2]) / vector[1] * vector[2],
#        center[1] - (y.limits[2]-center[2]) / vector[1] *vector[2])
# panel.points(x, y.limits, type="l", col="green")
            }
        }
    }
    
# plot points
    j <- 0
    for( i in include ){
        j <- j+1
        sel <- !flagFiltered & label==i
        panel.points(x[sel], y[sel], type="p", pch=pch[j], col=col[j], 
                    cex=cex[j])
    }
    
}

.clust.panel.splom <- function(x, y, frame, ...)
{
    cv <- current.viewport()$name
    dims <- as.numeric(strsplit(cv, ".", fixed=TRUE)[[1]][2:3])
    
    .clust.panel.xyplot(x=x, y=y, frame=frame, subset=c(dims[1],dims[2]), ...)
}

setMethod("splom", 
signature=signature(x="immunoClust", data="missing"),
definition=function(x, data, include=1:(x@K), ...)
{
    param <- attributes(x)$param
    
    y <- as.matrix(x@mu)
    colnames(y) <- param
    
    gp <- list(...)[["par.settings"]]
    
    ellipse.frame <- list(K=(x@K), P=length(param), sigma=(x@sigma), mu=(x@mu))
    
    splom(x=y, data=NULL, pscales=NULL, panel=.ellipse.panel.splom, 
        frame=ellipse.frame, gp=gp, include=include,...)
})

setMethod("splom", 
signature=signature(x="immunoClust",data="flowFrame"),definition=function(
x, data, include=1:(x@K), subset=1:length(attributes(x)$param), 
N=NULL,label=NULL, desc=NULL,...
) {
    params <- attributes(x)$param[subset]
    
    y <- .exprs(data, params)
    
    gp <- list(...)[["par.settings"]]
    
    par <- parameters(data)
    inc <- match(params, par@data[,'name'])
    range <- range(data)[inc]
    
    if( is.null(label) ){
        clust.frame <- list(K=(x@K), P=length(params), 
                            sigma=(x@sigma[,subset,subset]), 
                            mu=(x@mu[,subset]), range=range, label=x@label)
    }
    else {
        clust.frame <- list(K=(x@K), P=length(params), 
                            sigma=(x@sigma[,subset,subset]), 
                            mu=(x@mu[,subset]), range=range, label=label)
    }
    
    
    if( !is.null(N) && N < nrow(y) ) {
        y <- y[1:N,]
        clust.frame$label <- clust.frame$label[1:N]
    }
    
    varnames <- NULL
    
    name <- par@data[inc, 'name']
    if( is.null(desc) ) {
        if( is.null(attr(x, "desc")) )
        desc <- par@data[inc, 'desc']
        else
        desc <- attr(x,"desc")[inc]
    }
    else {
        desc <- desc[inc]
    }
    varnames <- paste(sep="\n", name, desc)
    
    splom(x=y, data=NULL, pscales=NULL, varnames=varnames, 
        panel=.clust.panel.splom, frame=clust.frame, gp=gp, 
        include=include, ...)
})

setMethod("splom", 
signature=signature(x="immunoClust",data="matrix"), definition=function(
x, data, include=1:(x@K), subset=1:length(attributes(x)$param), 
N=NULL,label=NULL, desc=NULL,...
){
    params <- attributes(x)$param[subset]
    
    gp <- list(...)[["par.settings"]]
        
    y <- data[, params]
    
    if( is.null(label) ) {
        label <- x@label
    }
    
    clust.frame <- list(K=(x@K), P=length(params), 
                    sigma=(x@sigma[,subset,subset]), mu=(x@mu[,subset]), 
                    range=NULL, label=label)
    
    
    dim(clust.frame$sigma) <- c(x@K, length(subset), length(subset))
    dim(clust.frame$mu) <- c(x@K, length(subset))
    
    if( !is.null(N) && N < nrow(y) ) {
        y <- y[1:N,]
        clust.frame$label <- clust.frame$label[1:N]
    }
    
    varnames <- NULL
    
    if( is.null(desc) ) {
        varnames <- params
    }
    else {
        varnames <- paste(sep="\n", params, desc)
    }
    
    splom(x=y, data=NULL, pscales=NULL, varnames=varnames, 
        panel=.clust.panel.splom, frame=clust.frame, gp=gp, 
        include=include, ...)
})

datSplom <- 
function(label, data, subset=1:ncol(data), include=1:nrow(data), ...) 
{
    param <- colnames(data)[subset]
    y <- data[, param]
    
    gp <- list(...)[["par.settings"]]
    
    
    clust.frame <- list(K=nrow(data), P=length(param), sigma=NULL, mu=NULL, 
                    range=NULL,label=label)
    
    splom(x=y, data=NULL, pscales=NULL, varnames=param, 
        panel=.clust.panel.splom, frame=clust.frame, gp=gp, 
        include=include, ...)
}

