
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
    if (small==2) xy.new[,2]=xy.new[,2]/ratio
    if (small==1) xy.new[,1]=xy.new[,1]/ratio
    xy.new + cbind(rep(loc[1],n), rep(loc[2],n))
}



setMethod("plot", signature(x="immunoClust", y="missing"),
function(x, data, subset=c(1,2), ellipse=TRUE, show.rm=FALSE, include=1:(x@K), 
main=NULL, col=include+1, pch=".", cex=0.6, 
col.rm=1, pch.rm=1, cex.rm=0.6, ecol=col, elty=1, 
npoints=501, add=FALSE, ...)
{
    
    if (!is.numeric(subset)) subset <- match(subset, x@parameters)
    
    if( !is.null(attr(data, "parameters")) ) {
        par <- parameters(data)
        inc <- match(x@parameters, par@data[,'name'])
        name <- par@data[inc, 'name']
        desc <- par@data[inc, 'desc']
        lab <- paste(sep=" / ", name, desc)[subset]
    }
    else {
        lab <- x@parameters[subset]
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
    
    if (!add) 
        plot(data, type="n", main=main, xlab=lab[1], ylab=lab[2], ...)  
    else 
        title(main)
    
    flagFiltered <- is.na(label)
# plot points with different colors/symbols corr. to cluster assignment
    col <- matrix(col, length(include))
    pch <- matrix(pch, length(include))
    cex <- matrix(cex, length(include))
    
    j <- 0
    
    for (i in include)  
        points(data[!flagFiltered & label==i,], 
                pch=pch[j <- j+1], col=col[j], cex=cex[j])  
    
# plot filtered points (from above or below)
    if (show.rm) 
        points(data[flagFiltered,], pch=pch.rm, col=col.rm, cex=cex.rm)
    
# plot ellipses
    if (ellipse) {
        ecol <- matrix(ecol, length(include))
        elty <- matrix(elty, length(include))
        
        cc <- qt(0.9, 5)
        
        j <- 0
        cc <- rep(cc, length.out=x@K)
        for (i in include) {
            eigenPair <- eigen(x@sigma[i,subset,subset])
            l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
            l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
            angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1]) 
            angle <- angle * 180/pi
            
            points(.ellipsePoints(a=l1[i], b=l2[i], alpha=angle, 
                                loc=x@mu[i,subset], n=npoints), 
                    type="l", lty=elty[j <- j+1], col=ecol[j])
        }  
    }
    
}
)

