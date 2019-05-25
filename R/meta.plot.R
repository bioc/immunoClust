## ploting

meta.plotClustersForScatter <- function(meta, include, filter=0)
{
    splom(meta$res.clusters, meta$dat.clusters$M, 
        include=.meta.ClustersForScatter(meta,include, filter=filter), 
        ellipse=TRUE)
}
meta.plotClusters <- function( meta, include=c() )
{
    if( !is.null(include) ) 
        splom(meta$res.clusters, meta$dat.clusters$M, include=include, 
        ellipse=TRUE)
    else
        splom(meta$res.clusters, meta$dat.clusters$M, ellipse=TRUE)
}

meta.plotScatter <- function(meta)
{
    if( ncol(meta$dat.scatter$M) == 2 )
    plot(meta$res.scatter, data=meta$dat.scatter$M, ellipse=TRUE)
    else
    splom(meta$res.scatter, meta$dat.scatter$M, ellipse=TRUE)
}

meta.plotGate <- function(meta, gating, pre, gates, pattern=c(), 
png.size=1024, filter=0)
{
    res <- meta$res.clusters 
    M <- meta$dat.clusters$M 
    name <- gsub("[/|\\]", "|", gating$desc)
    #cat(name, length(gating$par), "\n")
    
    if( "par" %in% attributes(gating)$names ) {
        
        gates[gating$par,] <- gating$gates

    }

    N <- length(pattern)

    if( N == 0 ) {
        png(filename=paste(sep="", pre, name, ".png"), width=png.size, 
            height=png.size)
        if( filter > 0 ) {
            inc <- c()
            for( c in gating$clusters ) {
                if( sum(res@label==c) > filter )
                inc <- c(inc,c)
            }
                
            if( length(inc) > 0 ) {
                print(splom(res, M, include=inc, ellipse=TRUE, gates=gates, 
                    ellipse.force=TRUE))
            }
        }
        else 
        if( !is.null(gating$clusters) ){
            print(splom(res, M, include=gating$clusters, ellipse=TRUE, 
                gates=gates, ellipse.force=TRUE))
        }
        dev.off()
    }

    childs <- gating$childs
    if( length(childs) > 0 ) {
        pre <- paste(sep="", pre, name, "_")
        for( g in seq_len(length(childs)) ) {
            if( N > 1 && childs[[g]]$desc == pattern[1] ) {
                meta.plotGate(meta, childs[[g]], pre, gates, 
                    pattern[2:N], png.size=png.size, filter=filter)
            }
            else
            if( N > 0 && childs[[g]]$desc == pattern[1] ) {
                meta.plotGate(meta, childs[[g]], pre, gates, 
                    png.size=png.size, filter=filter)
            }
            else 
            if( N==0 ){
                meta.plotGate(meta, childs[[g]], pre, gates, 
                    png.size=png.size, filter=filter)
            }
        }
    }
}

meta.plotGating <- function(meta, pre, pattern=c(), png.size=1024, filter=0)
{

    N <- length(pattern)
    if( N== 0 && !is.null(meta$res.scatter) ) {
        png(filename=paste(sep="", pre, "scatter.png"))
        print(meta.plotScatter(meta))
        dev.off()
    }

    childs <- meta$gating$childs
    gates <- matrix(NA, meta$dat.clusters$P, 3)

    if( length(childs) > 0 ) {
        for( g in seq_len(length(childs)) ) {
            if( N>1 && childs[[g]]$desc == pattern[1] ) {
                meta.plotGate(meta, childs[[g]], pre, gates, pattern[2:N], 
                png.size=png.size, filter=filter)
            }
            else
            if( N>0 && childs[[g]]$desc == pattern[1] ) {
                meta.plotGate(meta, childs[[g]], pre, gates, png.size=png.size, 
                    filter=filter)
            }
            else 
            if( N==0 ){
                meta.plotGate(meta, childs[[g]], pre, gates, png.size=png.size, 
                    filter=filter)
            }
        }
    }
}


meta.plotExpResult <- function(meta, exp, pattern=c(), png_pre, 
plot.ellipse=FALSE, plot.class=FALSE, png.size=1024, filter=0,
ellipse.quantile=0.95)
{

    gate <- .meta.getGate(meta$gating, pattern)

    if( is.null(gate) )
    return

    if( filter > 0 ) {
        inc <- c()
        for( c in gate$clusters ) {
            if( sum(meta$res.clusters@label==c) > filter )
            inc <- c(inc,c)
        }
    }
    else {
        inc <- gate$cluster
    }

    ##cat("inc", inc, "\n")

    cls <- c()
    cols <- c()
    j <- 1
    for( i in inc ) {
        j <- j+1
        cl <- which(meta$res.clusters@label==i)
        cls <- c(cls, cl )
        cols <- c(cols, rep(i+1,length(cl)))
    }

    for( i in seq_len(length(exp)) ) {
        res <- exp[[i]]
        k <- meta$dat.clusters$K[i]
        inc <- cls[ cls <= k ]
        col <- cols[ cls <= k ]
        cols <- cols[ cls > k ]
        cls <- cls[ cls > k ] - k
        if( length(inc) > 0 ) {
            #cat("plot", i, ";", inc, "\n")
            if( !plot.class ) {
                col <- inc+1 
            }

            fcs <- read.FCS(attr(res, "fcsName"))
            dat <- trans.ApplyToData(res, fcs)
            png(filename=paste(sep="", png_pre, attr(res, "expName"), ".png"), 
                width=png.size, height=png.size)
            w <- res@w[inc]
            m <- match(sort(w, decreasing=TRUE),w)
            print(splom(res, dat, include=inc[m], col=col[m], 
            ellipse=plot.ellipse, ellipse.quantile=ellipse.quantile, N=NULL))
            dev.off()
            
        }
    }

}
meta.plotExpClusters <- function(meta, exp, include=c(), png_pre, 
png.size=1024, plot.ellipse=FALSE, plot.class=FALSE, class.col=c(), filter=0,
ellipse.quantile=0.95, desc=NULL, N=NULL)
{

    if( is.null(include) )
    return

    if( filter > 0 ) {
        inc <- c()
        for( c in include ) {
            if( sum(meta$res.clusters@label==c) > filter )
            inc <- c(inc,c)
        }
    }
    else {
        inc <- include
    }

##  cat("inc", inc, "\n")

    cls <- c()
    cols <- c()
    j <- 1
    for( i in inc ) {
        cl <- which(meta$res.clusters@label==i)
        cls <- c(cls, cl )
        if( is.null(class.col) ) 
        cols <- c(cols, rep(i+1,length(cl)))
        else
        cols <- c(cols, rep(class.col[j], length(cl)))
        j <- j+1
    }
    
    for( i in seq_len(length(exp)) ) {
        res <- exp[[i]]
        k <- meta$dat.clusters$K[i]
        inc <- cls[ cls <= k ]
        col <- cols[ cls <= k ]
        cols <- cols[ cls > k ]
        cls <- cls[ cls > k ] - k
        if( length(inc) > 0 ) {
##cat("plot", i, ";", inc, "\n")
            if( !plot.class ) {
                col <- inc+1 
            }
            png_file <- paste(sep="", png_pre, attr(res, "expName"), ".png")
            if( !file.exists(png_file) ) {
                fcs <- read.FCS(attr(res, "fcsName"))
                dat <- trans.ApplyToData(res, fcs)
                N <- length(res@label)
                if( N < nrow(dat) ) {
                    dat <- dat[seq_len(N),]
                }
                png(filename=png_file, width=png.size, height=png.size)
                w <- res@w[inc]
                m <- match(sort(w, decreasing=TRUE),w)
                print(splom(res, dat, include=inc[m], N=N, desc=desc, 
                col=col[m], ellipse=plot.ellipse, 
                ellipse.quantile=ellipse.quantile))
                dev.off()
            }
            
        }
    }
}
###
