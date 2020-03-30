

.annotate.plotPop <- function(res, M, pop, pos, main="",
plot.childs=TRUE, plot.unclassified=FALSE, plot.parent=TRUE,
inc.childs=NULL, plot.subset=NULL, plot.ellipse=TRUE, plot.all=FALSE,
pscal=NULL,...)
{
    
    if( is.null(pop) ) {
        return
    }
    
    if( is.null(pscal) && !is.null(pop$pscales) ) 
        pscal = pop$pscales
    
    if( nchar(main) > 0 && nchar(pop$desc) > 0 )
    main <- paste(sep="_", main, pop$desc)
    else if (nchar(pop$desc) > 0 )
    main <- pop$desc
    
    if( length(pos) == 0 ) {
        if( !is.null(plot.subset) ) 
        subset <- plot.subset
        else
        if( !is.null(pop$plot.subset) ) 
        subset <- pop$plot.subset
        else
        subset=seq_len(length(attributes(res)$param))
        
        if( is.null(inc.childs) )
        inc.childs <- seq_len(length(pop$childs))
        
## collect clusters and colors
        if( plot.childs && !is.null(pop$plot.childs) && pop$plot.childs ) {
            if( length(pop$childs) > 0 ) {
                col <- c()
                clusters <- c()

                for( i in inc.childs ) {
                    clusters <- c(clusters, pop$childs[[i]]$clusters)
                    pop.color <- (i+1)
                    if( !is.null(pop$childs[[i]]$plot.color) ) {
                        pop.color <- pop$childs[[i]]$plot.color
                    }
                    col <- c(col,
                        rep(pop.color, length(pop$childs[[i]]$clusters)))
                }
                
                rest <- pop$clusters[ is.na( match(pop$clusters, clusters) ) ]
                if( plot.unclassified ) {
                    col <- c(rep("gray95", length(clusters)), rest + 1 )
                    clusters <- c(clusters, rest)
                }
                else {
                    clusters <- c(rest, clusters)
                    col <- c(rep("gray95", length(rest)), col)
                }
            }
            else {
                clusters <- pop$clusters
                if( !is.null(pop$plot.color) ) {
                    col <- rep(pop$plot.color, length(clusters))
                }
                else {
                    col <- clusters + 1  
                }
            }
        }
        else {
            clusters <- pop$clusters
            col <- clusters + 1
        }
        
        if( plot.parent && isTRUE(as.logical(pop$plot.parent)[1]) ) {
            parent.clusters <- c()
            if( is.logical(pop$plot.parent) )
                parent.clusters <- pop$parent.clusters
            else
                parent.clusters <- pop$plot.parent
                
            if( length(parent.clusters) > 0 ) {
                rest <- parent.clusters[is.na(match(parent.clusters, clusters))]
                clusters <- c(rest, clusters)
                col <- c(rep("gray95", length(rest)),col)
            }
        }

        mmain <- paste(sep=".", paste(pop$position, collapse="."), main)
        key <- NULL
        g <- length(pop$childs)
        if( g > 0 ) {
            key.pch <- rep(1,g)
            key.col <- seq_len(g)+1
            for( i in seq_len(g) ) if( !is.null(pop$childs[[i]]$plot.color) ) {
                key.col[i] <- pop$childs[[i]]$plot.color
            }
            #key.text <- sapply(pop$childs,function(x)x$desc)
            key.text <- vapply(pop$childs,function(x)x$desc, "")
            key <- list(title="Sub Level",
            cex.title=1,
            columns=1,
            points=list(pch=key.pch, col=key.col),
            text = list(key.text,adj=0, cex=1),
            adj=0,
            space="right",
            just="top")
        }
        
        if( length(subset) > 2 && length(clusters) > 0 ) {
            subpscal <- NULL
            
            if( !is.null(pscal) ) {
                subpscal <- list(seq_len((length(subset))))
                p <- 1
                for( par in subset ) {
                    subpscal[[p]] <- pscal[[par]]
                    p <- p+1
                }
            }
            
            if(plot.all) {
            print(splom(res, M, include=clusters, subset=subset, col=col,
                ellipse=plot.ellipse, main=mmain,
                mean=pop$M[subset], sigma=pop$S[subset, subset],
                pscales=subpscal, xlab=NULL, key=key, ... ))
            }
            else {
            return(splom(res, M, include=clusters, subset=subset, col=col,
                ellipse=plot.ellipse, main=mmain,
                mean=pop$M[subset], sigma=pop$S[subset, subset],
                pscales=subpscal, xlab=NULL, key=key, ... ))
            }
        }
        else
        if( length(clusters) > 0 ){
            plot(res, data=M, include=clusters, subset=subset, col=col, 
                ellipse=plot.ellipse, main=mmain, pscales=pscal,...)
            if( !is.null(key) ) {
            legend("topright", legend=key.text,pch=key.pch, col=key.col,
            title=key$title, title.adj=0, cex=1)
            }
        }
        
    } ## pos==null => plot
    
    if( length(pos) > 0 ) {
        p <- pos[1]
        
        if(length(pos) > 1 ) {
            pos <- pos[2:length(pos)]
        }
        else {
            pos <- NULL
        }
        
        .annotate.plotPop(res, M, pop$childs[[p]], pos, main=main,
        plot.childs=plot.childs,
        plot.unclassified=plot.unclassified,
        plot.parent=plot.parent,
        plot.ellipse=plot.ellipse,
        inc.childs=inc.childs,
        plot.subset=plot.subset,plot.all=plot.all,
        pscal=pscal,...)
    }
    else
    if( plot.all ) {
        #cat(paste(pop$position,collapse="."), pop$plot.endLevel, "\n")
        if( !is.null(pop$plot.endLevel) && pop$plot.endLevel ) {
            return(NULL)
        }
        
        for( i in seq_len(length(pop$childs)) ) {
        if( length(pop$childs[[i]]$childs) > 0 ||
            length(pop$childs[[i]]$clusters) > 1 ) {
            .annotate.plotPop( res, M, pop$childs[[i]], c(), main=main,
            plot.childs=plot.childs, plot.unclassified=plot.unclassified,
            plot.parent=plot.parent,
            plot.ellipse=plot.ellipse,
            plot.subset=plot.subset, plot.all=TRUE,
            pscal=pscal,...)
            }
        }
    }
}


setGeneric("plot")

plot.immunoMeta<-
function(x, pos=c(), main="", 
plot.childs=TRUE, plot.unclassified=FALSE,  plot.subset=c(),
inc.childs=c(), plot.ellipse=TRUE, plot.all=FALSE, ...)
{
    if( is.character(plot.subset) && plot.subset=="parent" ) {
        if(length(pos) > 0 )
        plot.subset <- prop(x, "plot.subset", pos[seq_len(length(pos)-1)])
        else
        plot.subset <- c()
    }
    .annotate.plotPop(x$res.clusters, x$dat.clusters$M, x$gating, pos,
    main=main,
    plot.childs=plot.childs,
    plot.unclassified=plot.unclassified,
    plot.subset=plot.subset, inc.childs=inc.childs,
    plot.ellipse=plot.ellipse,plot.all=plot.all, ...)
}
