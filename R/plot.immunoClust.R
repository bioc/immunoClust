
.ellipsePoints <- function(a,b, alpha = 0, loc = c(0,0), n = 501)
{
## Purpose: ellipse points,radially equispaced, given geometric par.s
## -------------------------------------------------------------------------
## Arguments: a, b : length of half axes in (x,y) direction
##            alpha: angle (in degrees) for rotation
##            loc  : center of ellipse
##            n    : number of points
## -------------------------------------------------------------------------
## Author: Martin Maechler, Date: 19 Mar 2002, 16:26
## modified by Kenneth to get rid of the precision problem met when 
## there s a large difference in the length of the two axes
    
    small <- 0
    if (a/b > 1000) {
        ratio <- a/b
        b <- a
        if (round(alpha)==0) small <- 2 else small <- 1
    }
    
    B <- min(a,b)
    A <- max(a,b)
## B <= A
    d2 <- (A-B)*(A+B)                   #= A^2 - B^2
    phi <- 2*pi*seq(0,1, len = n)
    sp <- sin(phi)
    cp <- cos(phi)
    r <- a*b / sqrt(B^2 + d2 * sp^2)
    xy <- r * cbind(cp, sp)
## xy are the ellipse points for alpha = 0 and loc = (0,0)
    al <- alpha * pi/180
    ca <- cos(al)
    sa <- sin(al)
    
    xy.new <- xy %*% rbind(c(ca, sa), c(-sa, ca))
    if (small==2) xy.new[,2]<-xy.new[,2]/ratio
    if (small==1) xy.new[,1]<-xy.new[,1]/ratio
    xy.new + cbind(rep(loc[1],n), rep(loc[2],n))
}

.plot_ellipses <- function(S, M, limits, cc=NULL, lty=3, col=NA, npoints=501)
{
    if( is.null(S) ) return()
    
    if( is.null(cc) )
        cc <- qt(0.95,5)
        
    if( is.null(limits) ) {
        eigenPair <- eigen(S)
        l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
        l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
        angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1])
        angle <- angle * 180/pi
        
        points(.ellipsePoints(a=l1, b=l2, alpha=angle,
                loc=M, n=npoints),
                type="l", lty=lty, col=col)
    }
    else {
        scale <- diag(1/limits)
        
        if( is.null(S) ) return ()
        
        S_scale <- scale %*% S %*% scale
        e <- NULL
        try( e <- eigen(S_scale) )
        if( is.null(e) ) {
            #cat(S, "\n", limits, "\n")
            return()
        }
        
        l1 <- sqrt(e$values[1]) * cc
        l2 <- sqrt(e$values[2]) * cc
        
        sa <- e$vectors[2,1]
        ca <- e$vectors[1,1]
        d2 <- (l1-l2)*(l1+l2)
        
        phi <- 2*pi*seq(0,1, len = npoints)
        sp <- sin(phi)
        cp <- cos(phi)
        
        r <- l1*l2 / sqrt(l2^2 + d2 * sp^2)
        xy <- r * cbind(cp, sp)
        xy <- xy %*% rbind(c(ca,sa), c(-sa,ca))
        #xy[,1] <- xy[,1] / scale
        xy <- xy %*% diag(limits)
        xy <- xy + cbind(rep(M[1],npoints), rep(M[2],npoints))
        
        polygon(xy[,1], xy[,2], border = col,  lty=lty, col=NA)
    }
}

setGeneric("plot")

setMethod("plot", signature(x="immunoClust", y="missing"),
function(x, data, subset=c(1,2), ellipse=TRUE, show.rm=FALSE, 
include=seq_len(x@K), main=NULL, col=include+1, pch=".", cex=0.6, 
col.rm=1, pch.rm=1, cex.rm=0.6, ecol=col, elty=1, 
npoints=501, add=FALSE, gates=NULL, pscales=NULL, ...)
{
    more.par <- list(...)
    
    if (!is.numeric(subset)) subset <- match(subset, x@parameters)
    
    if( !is.null(attr(x, "desc")) ) {
        desc <- attr(x, "desc")
        name <- gsub("\n", ": ", x@parameters)
        lab <- c()
        for( i in subset ) {
            if( !is.na(desc[i]) && desc[i] != "" )
            lab <- append(lab, paste(sep=": ", name[i], desc[i]) )
            else
            lab <- append(lab, name[i])
        }
        #lab <- paste(sep=": ", name, desc)[subset]
    }
    else
    if( !is.null(attr(data, "parameters")) ) {
        par <- parameters(data)
        inc <- match(x@parameters, par@data[,'name'])
        name <- par@data[inc, 'name']
        desc <- par@data[inc, 'desc']
        lab <- paste(sep=": ", name, desc)[subset]
    }
    else {
        lab <- gsub("\n", ": ", x@parameters[subset])
    }
    
    if (is(data, "flowFrame")) 
        data <- exprs(data)[,x@parameters] 
    else
    if (is(data, "matrix")) 
        (if (length(x@parameters)>0) data <- as.matrix(data[,x@parameters]))
    else
    if (is(data, "data.frame")) data <- as.matrix(data[,x@parameters])
    
    
    P <- ncol(data)
    data <- data[,subset]
    label <- x@label
    
    if( !is.null(more.par$N) && more.par$N < length(label) ) {
        label <- label[seq_len(more.par$N)]
    }
    
    limits <- NULL
    if (!add) {
        if( !is.null(pscales) ) {

            plot(data, type="n", main=main,
                xlab=lab[1], xlim=pscales[[ subset[1] ]]$limits,
                ylab=lab[2], ylim=pscales[[ subset[2] ]]$limits,
                xaxt="n", yaxt="n", ...)
            axis(1, at=pscales[[ subset[1] ]]$at,
                labels=pscales[[ subset[1] ]]$labels, ...)
            axis(2, at=pscales[[ subset[2] ]]$at,
                labels=pscales[[ subset[2] ]]$labels, ...)
            
            if( isTRUE(more.par$grid) ) {
                abline(h=pscales[[ subset[2] ]]$at, v=pscales[[ subset[1] ]]$at,
                    col="lightgray", lty=1)
            }
            cex.axis <- par("cex.axis")
            mtext(pscales[[ subset[1] ]]$unit, 1, line=3, adj=1, cex=cex.axis)
            mtext(pscales[[ subset[2] ]]$unit, 2, line=3, adj=1, cex=cex.axis)
            
            limits <- vapply(pscales, function(x)x$limits, numeric(2))
            
        }
        else {
            plot(data, type="n", main=main, xlab=lab[1], ylab=lab[2], ...)
        }
    }
    else {
        title(main)
    }
    
    flagFiltered <- is.na(label)
# plot points with different colors/symbols corr. to cluster assignment
    if( length(include) ) {
        col <- matrix(col, length(include))
        pch <- matrix(pch, length(include))
        cex <- matrix(cex, length(include))
    }
    else {
        col <- c()
        pch <- c()
        cex <- c()
    }
    
    j <- 0
    for (i in include) { 
        pts <- data[!flagFiltered & label==i,]
        j <- j+1
        if( is(pts, "matrix") ) {
        points(pts, pch=pch[j], col=col[j], cex=cex[j])
        }
        else {
        points(pts[1],pts[2], pch=pch[j], col=col[j], cex=cex[j])
        }
    }
# plot filtered points (from above or below)
    if (show.rm) 
        points(data[flagFiltered,], pch=pch.rm, col=col.rm, cex=cex.rm)
    
# plot ellipses
    if (ellipse) {
        ecol <- matrix(ecol, length(include))
        elty <- matrix(elty, length(include))
        
        cc <- qt(0.9, 5)
        j <- 0
        for (i in include) {
            j <- j+1
            .plot_ellipses(x@sigma[i,subset,subset], x@mu[i,subset],
                limits=limits[2,subset], cc=cc,
                lty <- elty[j], col=ecol[j], npoints=npoints)
        }
    }
    
    ## >> sollte mal verschwinden
    if( isTRUE(more.par$ellipsis.merged) ) {
        merged <- .clust.mergedClusters(x, include)
        
        .plot_ellipses(merged$sigma[subset,subset], merged$mu[subset],
            limits=limits[2,subset], cc=qt(0.95,5), lty=3, col="black",
            npoints=npoints)
    }
    
    if( !is.null(more.par$ellipses.mean) &&
        !is.null(more.par$ellipses.sigma) ) {
        
        if( is.null(more.par$ellipses.cc))
            cc <- rep(qt(0.98, 5), nrow(more.par$ellipses.mean))
        else
            cc <- more.par$ellipses.cc
            
        if( is.null(more.par$ellipses.col))
            col <- rep("black", nrow(more.par$ellipses.mean) )
        else
            col <- more.par$ellipses.col
            
        for( l in seq_len(nrow(more.par$ellipses.mean)) ) {
            
            .plot_ellipses(more.par$ellipses.sigma[l,subset,subset],
                more.par$ellipses.mean[l,subset], limits=limits[2,subset],
                cc=cc[l], lty=3, col=col[l], npoints=npoints)
        }
    }
    ## << weg damit

# plot gates    
    if( !is.null(gates) ) {
        x.limits <- c(min(data[!flagFiltered,1],-1),
                    max(data[!flagFiltered,1],10))
        y.limits <- c(min(data[!flagFiltered,2],-1),
                    max(data[!flagFiltered,2],10))
        thres <- gates[subset,]
        for( j in seq_len(length(thres[1,])) ) {   
            if( !is.na(thres[1,j]) ) {
                points(c(thres[1,j],thres[1,j]), y.limits,
                            type="l", col=j+1)
            }
        }
        for( j in seq_len(length(thres[2,])) ) {   
            if( !is.na(thres[2,j]) ) {
                points(x.limits, c(thres[2,j],thres[2,j]),
                            type="l", col=j+1)
            }
        }
    }
}
)

