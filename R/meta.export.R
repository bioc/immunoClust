####
### output tables
####
.cls_colors <- function()
{
c("black","red","green3","blue","cyan","magenta","yellow","gray")
}
###
## meta.numEvents: 
##        number of events in meta-clusters for cell-event samples
###
.meta.numGate <- function(meta, gating, pre, name, 
    out.linage=TRUE, out.all=TRUE, out.unclassified=TRUE)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    K <- meta$dat.clusters$K
    tbl <- matrix(NA, nrow=0, ncol=N+P)
    if( !is.null(gating) ) {
        
        if( is.null(name) ) {
            name <- gating$desc
        }
        else {
            name <- paste(sep="_", name, gating$desc)
        }

        clusters <- c()        
        if( out.linage && length(gating$clusters)>=0 ) {
            M  <- rep(0, P)
            inc <- c() 
            for( i in gating$clusters ) {
                j <- which(meta$res.clusters@label==i)
                M <- M + meta$res.clusters@mu[i,]*
                        sum(meta$dat.clusters$clsEvents[j])
                inc <- c(inc, j )
            }
            numEvents <- matrix(NA, nrow=1, ncol=N+P)
            if( length(inc) > 0 ) {
                M <- M / sum(meta$dat.clusters$clsEvents[inc])
                for( p in seq_len(P) ) {
                    numEvents[,N+p] <- M[p]
                }
            }
            
## get numEvents for all experiment
            k <- 0
            for( i in seq_len(N) ) {
                l <- k+1
                k <- k+K[i]
                ine <- inc[(inc>=l) & (inc<=k)]
                numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine])
            }
            
            if( out.all || length(gating$clusters) == 1 ) {
                tbl <- rbind(tbl, numEvents)
                if( out.unclassified && length(gating$clusters) == 1 )
                name <- paste(sep=".", name, gating$clusters, 
                            .cls_colors()[(gating$clusters%%8)+1])
                rownames(tbl) <- paste(sep=".", pre, name)
            }
            if( length(gating$clusters) == 1 )
                clusters <- c(clusters, gating$clusters)
            
            if( !is.null(gating$out_events) && nrow(gating$out_events) > 0 ) {
                numEvents <- matrix(NA, nrow=nrow(gating$out_events), 
                                        ncol=N+P)
                for( i in seq_len(N) ) {
                    for( j in seq_len(nrow(gating$out_events)) ) {
                        numEvents[j,i] <- gating$out_events[j,i]
                    }
                }
                rownames(numEvents) <- paste(sep="_", 
                                        paste(sep=".", pre, name), 
                                        rownames(gating$out_events))
                tbl <- rbind(tbl, numEvents)
            }
        }

        if( length(gating$childs) > 0 ) {
            c <- 0 
            for( i in seq_len(length(gating$childs)) ) 
                if( !is.null(gating$childs[[i]]) )
                    c <- c+1
            
            for( i in seq_len(length(gating$childs)) ) {
                clusters <- c(clusters, gating$childs[[i]]$clusters)
                if( is.null(pre) )
                pre_pos <- paste(sep="", i)
                else
                pre_pos <- paste(sep=".", pre, i)
                tbl <- rbind(tbl, .meta.numGate(meta, gating$childs[[i]], 
                                    pre_pos, name, 
                                    out.linage=(c>1||out.unclassified==FALSE), 
                                    out.all=out.all,
                                    out.unclassified=out.unclassified))
            }
            
        }
        clusters <- gating$clusters[ is.na(match(gating$clusters, clusters)) ]
        
        if( out.unclassified && length(clusters) > 0 ) {
            rn <- rownames(tbl)
            
            for( cls in clusters ) {
                inc <- which(meta$res.clusters@label==cls)
                numEvents <- matrix(NA, nrow=1, ncol=N+P)                
                k <- 0
                for( i in seq_len(N) ) {
                    l <- k+1
                    k <- k+K[i]
                    ine <- inc[(inc>=l) & (inc<=k)]
                    numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine])
                }
                for( p in seq_len(P) ) {
                    numEvents[,N+p] <- meta$res.clusters@mu[cls,p]
                }
                
                
                rn <- c(rn, paste(sep=".", pre, name, cls, 
                                    .cls_colors()[cls%%8+1]))
                tbl <- rbind(tbl, numEvents)
            }
            rownames(tbl) <- rn
        }
    }
    tbl
}
##    meta.numGate

meta.numEvents <- function(meta, out.all=TRUE, out.removed=FALSE, 
                    out.unclassified=TRUE)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    
    tbl <- matrix(NA, nrow=0, ncol=N+P)    
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
    
    
    if( out.removed ) {
        parEvents <- matrix(NA, nrow=2,ncol=N+P)
        parEvents[1,seq_len(N)] <- meta$dat.clusters$expEvents + 
                            meta$dat.clusters$removedEvents
        parEvents[2,seq_len(N)] <- meta$dat.clusters$removedEvents
        parNames <- c("measured", "removed")
        tbl <- rbind(tbl, parEvents)
        rownames(tbl) <- parNames
    }
    else
    if( out.all ) {
## all
        parEvents <- matrix(NA, nrow=1,ncol=N+P)    
        parEvents[1,seq_len(N)] <- meta$dat.clusters$expEvents
    
        parNames <- c("total")
    
        tbl <- rbind(tbl, parEvents)
        rownames(tbl) <- parNames
    }
    
    tbl <- rbind(tbl, .meta.numGate(meta, meta$gating, NULL, NULL,
                                    out.all=out.all, 
                                    out.unclassified=out.unclassified))
    
    colnames(tbl) <- c(meta$dat.clusters$expNames, parDesc)
    tbl
}
## meta.numEvents

###
## meta.numClusters: 
##        number of cell-clusters in meta-clusters for cell-event samples
###
.meta.clsGate <- function(meta, gating, name, out.linage=TRUE, out.all=TRUE)
{
    
##    cat("num clusters\n")
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    K <- meta$dat.clusters$K
    tbl <- matrix(NA, nrow=0, ncol=N)
    if( !is.null(gating) ) {
        
        if( is.null(name) ) {
            name <- gating$desc
        }
        else {
            name <- paste(sep="_", name, gating$desc)
        }

        clusters <- c()        
        if( out.linage && length(gating$clusters)>=0 ) {
            
            numClusters <- matrix(NA, nrow=1, ncol=N)
            
## get numEvents for all experiment
            inc <- c() 
            for( i in gating$clusters ) {
                j <- which(meta$res.clusters@label==i)
                inc <- c(inc, j )
            }
            k <- 0
            for( i in seq_len(N) ) {
                l <- k+1
                k <- k+K[i]
                numClusters[1,i] <- sum((inc>=l) & (inc<=k))
            }
            
            if( out.all || length(gating$clusters) == 1 ) {
                tbl <- rbind(tbl, numClusters)
                if( length(gating$clusters) == 1 )
                name <- paste(sep=".", name, gating$clusters, 
                            .cls_colors()[(gating$clusters%%8)+1])
                rownames(tbl) <- name
                if( length(gating$clusters) == 1 )
                    clusters <- c(clusters, gating$clusters)
            }
        }

        if( length(gating$childs) > 0 ) {
            c <- 0 
            for( i in seq_len(length(gating$childs)) ) 
            if( !is.null(gating$childs[[i]]) )
            c <- c+1
            
            for( i in seq_len(length(gating$childs)) ) {
                clusters <- c(clusters, gating$childs[[i]]$clusters)
                tbl <- rbind(tbl, .meta.clsGate(meta, gating$childs[[i]], name, 
                                        out.linage=c>1, out.all=out.all))
            }
            clusters <- gating$clusters[is.na(match(gating$clusters, clusters))]
            
        }
        if( length(clusters) > 0 ) {
            rn <- rownames(tbl)
            
            for( cls in clusters ) {
                inc <- which(meta$res.clusters@label==cls)
                numClusters <- matrix(NA, nrow=1, ncol=N)                
                k <- 0
                for( i in seq_len(N) ) {
                    l <- k+1
                    k <- k+K[i]
                    numClusters[1,i] <- sum((inc>=l) & (inc<=k))
                }
                
                rn <- c(rn, paste(sep=".", name, cls, .cls_colors()[cls%%8+1]))
                tbl <- rbind(tbl, numClusters)
            }
            rownames(tbl) <- rn
        }
    }
    tbl
}
##    meta.clsGate

meta.numClusters <- function(meta, out.all=TRUE)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    
    tbl <- matrix(NA, nrow=0, ncol=N)    
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
    
    
## all
    parEvents <- matrix(NA, nrow=1,ncol=N)
    
    parEvents[1,seq_len(N)] <- meta$dat.clusters$K
    
    parNames <- c("total")
    
    if( out.all ) {
        tbl <- rbind(tbl, parEvents)
        rownames(tbl) <- parNames
    }
    
    tbl <- rbind(tbl, .meta.clsGate(meta, meta$gating, NULL, out.all=out.all))
    
    
    colnames(tbl) <- c(meta$dat.clusters$expNames)
    tbl
}
## meta.numClusters


###
## meta.relEvents: 
##        relative event frequencies (per total events)
###
meta.relEvents <- function(meta, out.all=TRUE, out.removed=FALSE,
                    out.unclassified=TRUE)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    
    tbl <- matrix(NA, nrow=0, ncol=N+P)    
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
    
    
## all
    parEvents <- matrix(NA, nrow=1,ncol=N+P)
    
    
    if( out.removed ) {
        parEvents[1,seq_len(N)] <- meta$dat.clusters$removedEvents
        rownames(parEvents) <- "removed"
        tbl <- rbind(tbl, parEvents)
    }
    
    if( out.all ) {
        parEvents[1,seq_len(N)] <- meta$dat.clusters$expEvents
        rownames(parEvents) <- "total"
        tbl <- rbind(tbl, parEvents)
    }
    
    tbl <- rbind(tbl, .meta.numGate(meta, meta$gating, NULL, NULL, 
                                    out.all=out.all,
                                    out.unclassified=out.unclassified))
    
    parEvents[1,seq_len(N)] <- meta$dat.clusters$expEvents +
                            meta$dat.clusters$removedEvents
    rownames(parEvents) <- "measured"
    if( nrow(tbl) > 0 ) {
        for( i in seq_len(nrow(tbl)) ) {
            for( j in seq_len(N) ) {
                tbl[i,j] <- tbl[i,j] / parEvents[1,j]
            }
        }
    }
    
    colnames(tbl) <- c(meta$dat.clusters$expNames, parDesc)
    tbl
}
## meta.relEvents

###
## meta.parMFI: 
###
.meta.parGate <- 
function(meta, gating, par, pre, name, out.linage=TRUE, out.all=TRUE)
{    
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    K <- meta$dat.clusters$K
    tbl <- matrix(NA, nrow=0, ncol=N+P)
    if( !is.null(gating) ) {
        
        if( is.null(name) ) {
            name <- gating$desc
        }
        else {
            name <- paste(sep="_", name, gating$desc)
        }

        clusters <- c() 
        if( out.linage && length(gating$clusters)>=0 ) {
            M  <- rep(0, P)
            inc <- c() 
            for( i in gating$clusters ) {
                j <- which(meta$res.clusters@label==i)
                M <- M + meta$res.clusters@mu[i,]*
                        sum(meta$dat.clusters$clsEvents[j])
                
                inc <- c(inc, j )
            }
            numEvents <- matrix(NA, nrow=1, ncol=N+P)
            
            if( length(inc) > 0 ) {
                M <- M / sum(meta$dat.clusters$clsEvents[inc])
            
                for( p in seq_len(P) ) {
                    numEvents[,N+p] <- M[p]
                }
            
            
## get numEvents for all experiment
                k <- 0
                for( i in seq_len(N) ) {
                l <- k+1
                k <- k+K[i]
                ine <- inc[(inc>=l) & (inc<=k)]
## weighted MFI center
                numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine]*
                                    meta$dat.clusters$M[ine,par])/
                                    sum(meta$dat.clusters$clsEvents[ine])
                }
            }
            
            if( out.all || length(gating$clusters) == 1 ) {
                tbl <- rbind(tbl, numEvents)
                if( length(gating$clusters) == 1 )
                name <- paste(sep=".", name, gating$clusters, 
                            .cls_colors()[(gating$clusters%%8)+1])
                rownames(tbl) <- paste(sep=".", pre, name)
                if( length(gating$clusters) == 1 ) 
                    clusters <- c(clusters, gating$clusters)
            }
        }
        

        if( length(gating$childs) > 0 ) {
            c <- 0 
            for( i in seq_len(length(gating$childs)) ) 
                if( !is.null(gating$childs[[i]]) )
                    c <- c+1
            
            for( i in seq_len(length(gating$childs)) ) {
                clusters <- c(clusters, gating$childs[[i]]$clusters)
                if( is.null(pre) )
                pre_pos <- paste(sep="", i)
                else 
                pre_pos <- paste(sep=".", pre, i)
                
                tbl <- rbind(tbl, .meta.parGate(meta, gating$childs[[i]], 
                                        par, pre_pos, name, 
                                        out.linage=c>1, out.all=out.all))
            }
            
        }
        clusters <- gating$clusters[ is.na(match(gating$clusters, clusters)) ]
        
        if( length(clusters) > 0 ) {
            rn <- rownames(tbl)

            for( cls in clusters ) {
                inc <- which(meta$res.clusters@label==cls)
                numEvents <- matrix(NA, nrow=1, ncol=N+P) 
                k <- 0
                for( i in seq_len(N) ) {
                    l <- k+1
                    k <- k+K[i]
                    ine <- inc[(inc>=l) & (inc<=k)]
## weightes MFI center
                    numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine]*
                                        meta$dat.clusters$M[ine,par])/
                                        sum(meta$dat.clusters$clsEvents[ine])
                }
                for( p in seq_len(P) ) {
                    numEvents[,N+p] <- meta$res.clusters@mu[cls,p]
                }
                
                cls_name <- paste(sep=".", pre, name, cls, 
                                    .cls_colors()[cls%%8+1])
                rn <- c(rn, cls_name)
                tbl <- rbind(tbl, numEvents)
            }
            rownames(tbl) <- rn
        }
    }
    tbl
}
##
#    meta.parGate
##

meta.parMFI <- function(meta, par, out.all=TRUE)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    K <- meta$dat.clusters$K
    
    tbl <- matrix(NA, nrow=0, ncol=N+P)    
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
    
    
## all
    numEvents <- matrix(NA, nrow=1,ncol=N+P)
    k <- 0
    for( i in seq_len(N) ) {
        l <- k+1
        k <- k+K[i]
        
        numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[l:k]*
                            meta$dat.clusters$M[l:k,par])/
                            sum(meta$dat.clusters$clsEvents[l:k])
    }
    
    if( out.all ) {
        tbl <- rbind(tbl, numEvents)
        rownames(tbl) <- c("total")
    }
    
    tbl <- rbind(tbl, .meta.parGate(meta, meta$gating, par, NULL, NULL, 
                                out.all=out.all))
    
    colnames(tbl) <- c(meta$dat.clusters$expNames, parDesc)
    tbl
}
## meta.parMFI

###
## meta.freqTable:
###
.meta.freqGate <- function(meta, gating, totEvents, name, out.linage=TRUE)
{
    
##    cat("gate frequencies:", nrow(totEvents), "\n")
    N <- meta$dat.clusters$N
    K <- meta$dat.clusters$K
    P <- meta$dat.clusters$P
    
    tbl <- matrix(NA, nrow=0, ncol=N+P)
    parEvents <- totEvents
    if( !is.null(gating) ) {
        
        if( is.null(name) ) {
            name <- gating$desc
        }
        else {
            name <- paste(sep="_", name, gating$desc)
        }
        
        if( out.linage && length(gating$clusters)>0 ) {
            
            inc <- c()
            M  <- rep(0, P)
            for( i in gating$clusters ) {
                j <- which(meta$res.clusters@label==i)
                M <- M + meta$res.clusters@mu[i,]*
                        sum(meta$dat.clusters$clsEvents[j])
                inc <- c(inc, j )
            }
            
            M <- M / sum(meta$dat.clusters$clsEvents[inc])
            
            numEvents <- matrix(NA, nrow=nrow(totEvents), ncol=N+P)
            if( length(gating$cluster) == 1 ){
                name <- paste(sep=".", name, gating$clusters, 
                            .cls_colors()[(gating$clusters%%8)+1])
                
                for( p in seq_len(P) ) {
                    numEvents[,N+p] <- meta$res.clusters@mu[gating$clusters,p]
                }
            }
            else {
## weighted sum
                for( p in seq_len(P) ) {
                    numEvents[,N+p] <- M[p]
                }
                
            }
            
## get frequencies for all experiment
            k <- 0
            for( i in seq_len(N) ) {
                l <- k+1
                k <- k+K[i]
                ine <- inc[(inc>=l) & (inc<=k)]
                for( j in seq_len(nrow(totEvents)) ) {
                    numEvents[j,i] <- sum(meta$dat.clusters$clsEvents[ine]) / 
                                        totEvents[j,i]
                }
            }
            
            tbl <- rbind(tbl, numEvents)
            rownames(tbl) <- paste(sep="", name, "./.", rownames(totEvents))
## num events             
            numEvents <- matrix(NA, nrow=1, ncol=N)
            k <- 0
            for( i in seq_len(N) ) {
                l <- k+1
                k <- k+K[i]
                ine <- inc[(inc>=l) & (inc<=k)]
                numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine])
            }
            
            rownames(numEvents) <- name
            parEvents <- rbind(totEvents, numEvents)
            
        }
        
        if( length(gating$childs) > 0 ) {
            
            c <- 0 
            for( i in seq_len(length(gating$childs)) ) 
            if( !is.null(gating$childs[[i]]) )
            c <- c+1
            
            for( i in seq_len(length(gating$childs)) ) 
            tbl <- rbind(tbl, .meta.freqGate(meta, gating$childs[[i]], 
                                            parEvents, name, out.linage=c>1))
            
        }
        else
        if( length(gating$clusters) > 1 ) {
            rn <- rownames(tbl)
            
            for( cls in gating$clusters ) {
                inc <- which(meta$res.clusters@label==cls)
                numEvents <- matrix(NA, nrow=1, ncol=N)                
                freqs <- matrix(NA, nrow=nrow(parEvents), ncol=N+P) 
                
                k <- 0
                for( i in seq_len(N) ) {
                    l <- k+1
                    k <- k+K[i]
                    ine <- inc[(inc>=l) & (inc<=k)]
                    numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine])
                }
                
                for( j in seq_len(nrow(parEvents)) ) {                
                    freqs[j,seq_len(N)] <- numEvents[1,] / parEvents[j,]
                    rn <- c(rn, paste(sep=".", name, cls, 
                                    .cls_colors()[cls%%8+1], 
                                    "/", rownames(parEvents)[j]))
                }
                
                for( p in seq_len(P) ) {
                    freqs[,N+p] <- meta$res.clusters@mu[cls,p]
                }
                
                tbl <- rbind(tbl, freqs)
            }
            rownames(tbl) <- rn
        }
    }
    tbl
}
##    meta.freqGate

meta.freqTable <- function(meta)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    
    tbl <- matrix(NA, nrow=0, ncol=N+P)    
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
## all
    parEvents <- matrix(meta$dat.clusters$expEvents, nrow=1,ncol=N)
    parNames <- c("total")
    rownames(parEvents) <- parNames
    
    tbl <- rbind(tbl, .meta.freqGate(meta, meta$gating, parEvents, NULL))
    
    colnames(tbl) <- c(meta$dat.clusters$expNames, parDesc)
    tbl
}
### meta.freqTable

### meta.relParent
.meta.numEvents <- function(res, dat, clusters)
{
    N <- dat$N
    P <- dat$P
    K <- dat$K
    
    M  <- rep(0, P)
    inc <- c()
    for( i in clusters ) {
        j <- which(res@label==i)
        #if( i < 1 || i > nrow(res@mu))
        #cat("out of bound", i, ":", clusters, "\n")
        M <- M + res@mu[i,] * sum(dat$clsEvents[j])
        inc <- c(inc, j )
    }
    numEvents <- matrix(NA, nrow=1, ncol=N+P)
    if( length(inc) > 0 ) {
        M <- M / sum(dat$clsEvents[inc])
        
        for( p in seq_len(P) ) {
            numEvents[,N+p] <- M[p]
        }
    }
    
    
    ## get numEvents for all experiment
    k <- 0
    for( i in seq_len(N) ) {
        l <- k+1
        k <- k+K[i]
        ine <- inc[(inc>=l) & (inc<=k)]
        numEvents[1,i] <- sum(dat$clsEvents[ine])
    }
    
    numEvents
}
.meta.relParent <- function(meta, gating, pre, name,
    out.linage=TRUE, out.all=TRUE, out.unclassified=TRUE)
{
    #cat(gating$desc, ":", gating$clusters, "\n")
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    K <- meta$dat.clusters$K
    tbl <- matrix(NA, nrow=0, ncol=N+P)
    if( !is.null(gating) ) {
        
        if( is.null(name) ) {
            name <- gating$desc
        }
        else {
            name <- paste(sep="_", name, gating$desc)
        }
        if( !is.null(gating$parent.clusters) ||
            length(gating$parent.position) > 1) {
            parEvents <- .meta.numEvents(meta$res.clusters, meta$dat.clusters,
            gating$parent.clusters)
            parDesc <- gating$parent.desc
        }
        else {
            parEvents <- matrix(NA, nrow=1, ncol=N+P)
            parEvents[1,seq_len(N)] <- (meta$dat.clusters$expEvents +
                                meta$dat.clusters$removedEvents)
            parDesc <- "measured"
        }
        
        ##cat( paste(collapse=".", gating$position), parDesc, "\n")
        clusters <- c()
        classified <- .annotate.retrieveClassified(gating,c())
        if( out.linage && length(gating$clusters)>=0 ) {
            
            numEvents <- .meta.numEvents(meta$res.clusters, meta$dat.clusters,
            gating$clusters)
            
            numEvents[1,seq_len(N)] <- numEvents[1,seq_len(N)] / 
                                        parEvents[1,seq_len(N)]
            
            
            if( out.all || length(gating$clusters) == 1 ) {
                tbl <- rbind(tbl, numEvents)
                
                if( length(gating$clusters) == 1 && out.unclassified
                && match(gating$clusters, classified,nomatch=0)==0 ){
                    cls_name <- paste(sep=".", name, gating$clusters,
                    .cls_colors()[(gating$clusters%%8)+1])
                }
                else {
                    cls_name <- name
                }
                
                if( !is.null(parEvents) && nchar(parDesc) > 0 ) {
                    cls_name <- sprintf("%s [/%s]", cls_name, parDesc)
                }
                rownames(tbl) <- paste(sep=".", pre, cls_name)
                if( length(gating$clusters) == 1 )
                clusters <- c(clusters, gating$clusters)
                
                
                if( !is.null(gating$out_events) &&
                    nrow(gating$out_events) > 0 ) {
                    clsEvents <- .meta.numEvents(meta$res.clusters,
                                    meta$dat.clusters, gating$clusters)
                    
                    numEvents <- matrix(NA, nrow=nrow(gating$out_events),
                    ncol=N+P)
                    for( i in seq_len(N) ) {
                        for( j in seq_len(nrow(gating$out_events)) ) {
                            numEvents[j,i] <- gating$out_events[j,i] / 
                                                clsEvents[1,i]
                        }
                    }
                    rn <- paste(sep="_", paste(sep=".", pre, name),
                    rownames(gating$out_events))
                    rn <- sprintf("%s [/%s]", rn, gating$desc)
                    rownames(numEvents) <- rn
                    tbl <- rbind(tbl, numEvents)
                }
                
            }
        }
        
        
        if( length(gating$childs) > 0 ) {
            c <- 0
            for( i in seq_len(length(gating$childs)) )
            if( !is.null(gating$childs[[i]]) )
            c <- c+1
            
            for( i in seq_len(length(gating$childs)) ) {
                clusters <- c(clusters, gating$childs[[i]]$clusters)
                if( is.null(pre) )
                pre_pos <- paste(sep="", i)
                else
                pre_pos <- paste(sep=".", pre, i)
                tbl <- rbind(tbl,
                            .meta.relParent(meta, gating$childs[[i]],
                                pre_pos, name,
                                out.linage=(c>1 || out.unclassified==FALSE),
                                out.all=out.all,
                                out.unclassified=out.unclassified)
                            )
            }
            
        }
        if( out.unclassified ) {
            clusters <- gating$clusters[is.na(match(gating$clusters, clusters))]
            
            if( length(clusters) > 0 ) {
                rn <- rownames(tbl)
                if( !is.null(gating$unclassified.parent.clusters) ) {
                    parEvents <- .meta.numEvents(meta$res.clusters,
                                    meta$dat.clusters,
                                    gating$unclassified.parent.clusters)
                    parDesc <- gating$unclassified.parent.desc
                }
                
                if( is.null(parEvents) ) {
                    parEvents <- .meta.numEvents(meta$res.clusters,
                                    meta$dat.clusters,
                                    gating$clusters)
                    parDesc <- gating$desc
                }
                #else {
                #    parDesc <- gating$parent.desc
                #}
                
                for( cls in clusters ) {
                    numEvents <- .meta.numEvents(meta$res.clusters,
                                    meta$dat.clusters, cls)
                    cls_name <- paste(sep=".", pre, name, cls,
                                        .cls_colors()[cls%%8+1])
                    if( !is.null(parEvents) ) {
                        numEvents[1,seq_len(N)] <- numEvents[1,seq_len(N)] / 
                                                    parEvents[1,seq_len(N)]
                        cls_name <- sprintf("%s, [/%s]", cls_name, parDesc)
                    }
                    rn <- c(rn, cls_name)
                    tbl <- rbind(tbl, numEvents)
                }
                rownames(tbl) <- rn
            }
        }
    }
    tbl
}

meta.relParent <-
function(meta, out.all=TRUE, out.removed=FALSE, out.unclassified=TRUE)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    K <- meta$dat.clusters$K
    
    tbl <- matrix(NA, nrow=0, ncol=N+P)
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
    
    
    if( out.removed ) {
        parEvents <- matrix(NA, nrow=2,ncol=N+P)
        parEvents[1,seq_len(N)] <- meta$dat.clusters$expEvents +
                                meta$dat.clusters$removedEvents
        parEvents[2,seq_len(N)] <- meta$dat.clusters$removedEvents / 
                                    parEvents[1,seq_len(N)]
        parNames <- c("measured", "removed [/measured]")
        tbl <- rbind(tbl, parEvents)
        rownames(tbl) <- parNames
        
    }
    else
    if( out.all ) {
        parEvents <- matrix(NA, nrow=1,ncol=N+P)
        parEvents[1,seq_len(N)] <- meta$dat.clusters$expEvents
        parNames <- c("total")
        
        tbl <- rbind(tbl, parEvents)
        rownames(tbl) <- parNames
    }
    
    tbl <- rbind(tbl,
                .meta.relParent(meta, meta$gating, NULL, NULL,
                out.all=out.all,
                out.unclassified=out.unclassified)
                )
    
    
    colnames(tbl) <- c(meta$dat.clusters$expNames, parDesc)
    tbl
}
## meta.relParent

################################################################################
