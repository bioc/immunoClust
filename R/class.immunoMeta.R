

immunoMeta <- function(res, dat, gating) {
    if( missing(dat) || is.null(dat) ) {
        clsEvents <- rep(0, res@K)
        for( k in seq_len(res@K) ) {
            clsEvents[k] <- sum(!is.na(res@label) & res@label==k)
        }
        M <- res@mu
        desc <- attr(res, "desc")
        if( !is.null(desc) ) {
            colnames(M) <- paste(sep="\n", res@parameters, desc)
            res@parameters <- colnames(M)
        }
        else {
            desc <- res@parameters
            colnames(M) <- attr(res, "parameters")
        }

        
        dat <- list("P"=res@P, "N"=1, "K"=res@K, 
                    "W"=res@w, "M"=M, "S"=res@sigma,
                    "expNames"=res@expName, 
                    "expEvents"=sum(!is.na(res@label)),
                    "removedEvents"=sum(is.na(res@label)), 
                    "clsEvents"=clsEvents, "desc"=desc)
        res@label <- seq_len(res@K)
        
    }
    if( missing(gating) || is.null(gating) ) {
        gating <- list("clusters"=seq_len(res@K), "childs"=c(), 
                        "desc"="all", "partition"=TRUE)
    }
    #names(res@w) <- sprintf("cls-%d", seq_len(res@K))
    structure(list("res.clusters"=res, "dat.clusters"=dat, "gating"=gating), 
        class = "immunoMeta")
    
}

setClass("immunoMeta")

setGeneric("summary")
summary.immunoMeta <-
function(object, ...){
    
    .level.summary <- function(pop, idx, depth=0)
    {
        if( is.null(pop) ) {
            ##cat("\tempty\n")
            return (c())
        }
        
        cls <- c()
        if( length(pop$childs) > 0 ) {
            
            for( i in seq_len(length(pop$childs)) ) {
                if( !is.null(pop$childs[[i]]) ) {
                    cls <- c(cls, pop$childs[[i]]$clusters)
                }
            }
            
            cls <- unique(cls)
        }
        plot.col <- idx+1
        if( !is.null(pop$plot.color) ) {
            plot.col <- pop$plot.color
        }
        if( is.numeric(plot.col)) {
            if( (plot.col <- plot.col%%8) == 0 )
            plot.col <- 8
            
            plot.col <- .cls_colors()[plot.col]
            
        }
        
        cat(paste(pop$position, collapse="."),
            descFull(object,pop$position), "\n")
        cat(sprintf("\t%d sub-levels, %d clusters, %s sub-classified\n",
            length(pop$childs),length(pop$clusters), length(cls)))

        cat("\tparent: ", paste(pop$parent.position, collapse="."), "\n")
        cat("\tplot:   ", plot.col, paste(pop$plot.subset,collapse=","), "\n")
        
        
        if( length(pop$childs) > 0 ) {
            for( i in seq_len(length(pop$childs)) ) {
                .level.summary( pop$childs[[i]], i, depth=depth+1)
            }
        }
        
    }
    
    cat("** Experiment Information ** \n")
    cat("Samples:     ", object$dat.clusters$N,"\n")
    cat("\tclusters    events  name\n")
    for( n in seq_len(object$dat.clusters$N)) {
        cat(sprintf("\t%8d  %8d  %s\n", object$dat.clusters$K[n],
            round(object$dat.clusters$expEvents[n]),
            object$dat.clusters$expNames[n]))
    }
    cat("Parameters:  ", npar(object), "\n")
    for( p in seq_len(npar(object)) ) {
        cat("\t", gsub("\n", ": ", parameters(object)[p]),"\n")
    }
    
    cat("\n** Clustering Summary ** \n")
    cat("Number of observations:", nobs(object),"\n")
    cat("Number of parameters:  ", npar(object),"\n")
    cat("ICL bias:              ", formatC(attr(object$res.clusters,"bias"),
                                            format="f",digits=2),"\n")
    cat("Number of clusters:    ", ncls(object),"\n")
    
    cat("\n** Annotation Summary ** \n")
    .level.summary(object$gating, 0)
}
setMethod("summary", signature(object="immunoMeta"),
function(object)
{
    summary.immunoMeta(object)
})

setGeneric("show")
show.immunoMeta <- function(object)
{
    cat("Object of class 'immunoMeta'\n")
    cat("This object has the following slots:\n")
    cat("$res.clusters 'immunoClust' object of the meta-clustering result\n" )
    cat("$dat.clusters 'list' object containing information of the clustered 
        cell-clusters\n")
    cat("$gating       'list' object containing the meta-clusters annotation 
        tree\n")
}

setMethod("show", signature(object="immunoMeta"),
function(object)
{
    show.immunoMeta(object)
})

