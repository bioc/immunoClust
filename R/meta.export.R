####
### output tables
####

###
## meta.numEvents: 
##        number of events in meta-clusters for cell-event samples
###
.meta.numGate <- function(meta, gating, name, out.linage=TRUE, out.all=TRUE)
{
    
#    cat("num events\n")
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
        
        if( out.linage && length(gating$clusters)>0 ) {
            M  <- rep(0, P)
            inc <- c() 
            for( i in gating$clusters ) {
                j <- which(meta$res.clusters@label==i)
                M <- M + meta$res.clusters@mu[i,]*
                        sum(meta$dat.clusters$clsEvents[j])
                
                inc <- c(inc, j )
            }
            M <- M / sum(meta$dat.clusters$clsEvents[inc])
            
            numEvents <- matrix(NA, nrow=1, ncol=N+P)
            for( p in 1:P ) {
                numEvents[,N+p] <- M[p]
            }
            
## get numEvents for all experiment
            k <- 0
            for( i in 1:N ) {
                l <- k+1
                k <- k+K[i]
                ine <- inc[(inc>=l) & (inc<=k)]
                numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine])
            }
            
            if( out.all || length(gating$clusters) == 1 ) {
                tbl <- rbind(tbl, numEvents)
                if( length(gating$clusters) == 1 )
                name <- paste(sep=".", name, gating$clusters, 
                            palette()[(gating$clusters%%8)+1])
                rownames(tbl) <- name
            }
        }
        
        if( !is.null(gating$childs) ) {
            c <- 0 
            for( i in 1:length(gating$childs) ) 
            if( !is.null(gating$childs[[i]]) )
            c <- c+1
            
            for( i in 1:length(gating$childs) ) 
            tbl <- rbind(tbl, .meta.numGate(meta, gating$childs[[i]], name, 
                                        out.linage=c>1, out.all=out.all))
            
        }
        else
        if( length(gating$clusters) > 1 ) {
            rn <- rownames(tbl)
            
            for( cls in gating$clusters ) {
                inc <- which(meta$res.clusters@label==cls)
                numEvents <- matrix(NA, nrow=1, ncol=N+P)                
                k <- 0
                for( i in 1:N ) {
                    l <- k+1
                    k <- k+K[i]
                    ine <- inc[(inc>=l) & (inc<=k)]
                    numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine])
                }
                for( p in 1:P ) {
                    numEvents[,N+p] <- meta$res.clusters@mu[cls,p]
                }
                
                
                rn <- c(rn, paste(sep=".", name, cls, palette()[cls%%8+1]))
                tbl <- rbind(tbl, numEvents)
            }
            rownames(tbl) <- rn
        }
    }
    tbl
}
##    meta.numGate

meta.numEvents <- function(meta, out.all=TRUE)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    
    tbl <- matrix(NA, nrow=0, ncol=N+P)    
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
    
    
## all
    parEvents <- matrix(NA, nrow=1,ncol=N+P)
    
    parEvents[1,1:N] <- meta$dat.clusters$expEvents
    
    parNames <- c("total")
    
    if( out.all ) {
        tbl <- rbind(tbl, parEvents)
        rownames(tbl) <- parNames
    }
    
    tbl <- rbind(tbl, .meta.numGate(meta, meta$gating, NULL, out.all=out.all))
    
    
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
        
        if( out.linage && length(gating$clusters)>0 ) {
            
            numClusters <- matrix(NA, nrow=1, ncol=N)
            
## get numEvents for all experiment
            inc <- c() 
            for( i in gating$clusters ) {
                j <- which(meta$res.clusters@label==i)
                inc <- c(inc, j )
            }
            k <- 0
            for( i in 1:N ) {
                l <- k+1
                k <- k+K[i]
                numClusters[1,i] <- sum((inc>=l) & (inc<=k))
            }
            
            if( out.all || length(gating$clusters) == 1 ) {
                tbl <- rbind(tbl, numClusters)
                if( length(gating$clusters) == 1 )
                name <- paste(sep=".", name, gating$clusters, 
                            palette()[(gating$clusters%%8)+1])
                rownames(tbl) <- name
            }
        }
        
        if( !is.null(gating$childs) ) {
            c <- 0 
            for( i in 1:length(gating$childs) ) 
            if( !is.null(gating$childs[[i]]) )
            c <- c+1
            
            for( i in 1:length(gating$childs) ) 
            tbl <- rbind(tbl, .meta.clsGate(meta, gating$childs[[i]], name, 
                                        out.linage=c>1, out.all=out.all))
            
        }
        else
        if( length(gating$cluster) > 1 ) {
            rn <- rownames(tbl)
            
            for( cls in gating$clusters ) {
                inc <- which(meta$res.clusters@label==cls)
                numClusters <- matrix(NA, nrow=1, ncol=N)                
                k <- 0
                for( i in 1:N ) {
                    l <- k+1
                    k <- k+K[i]
                    numClusters[1,i] <- sum((inc>=l) & (inc<=k))
                }
                
                rn <- c(rn, paste(sep=".", name, cls, palette()[cls%%8+1]))
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
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parmeters)
    
    
## all
    parEvents <- matrix(NA, nrow=1,ncol=N)
    
    parEvents[1,1:N] <- meta$dat.clusters$K
    
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
meta.relEvents <- function(meta, out.all=TRUE)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    
    tbl <- matrix(NA, nrow=0, ncol=N+P)    
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
    
    
## all
    parEvents <- matrix(NA, nrow=1,ncol=N+P)
    
    parEvents[1,1:N] <- meta$dat.clusters$expEvents
    
    parNames <- c("total")
    
    if( out.all ) {
        tbl <- rbind(tbl, parEvents)
        rownames(tbl) <- parNames
    }
    
    tbl <- rbind(tbl, .meta.numGate(meta, meta$gating, NULL, out.all=out.all))
    
    for( i in 1:nrow(tbl) ) {
        for( j in 1:N ) {
            tbl[i,j] <- tbl[i,j] / parEvents[1,j]
        }
    }
    
    colnames(tbl) <- c(meta$dat.clusters$expNames, parDesc)
    tbl
}
## meta.relEvents

meta.relEvents2 <- function(meta, major=1:5, out.all=TRUE)
{
    N <- meta$dat.clusters$N
    P <- meta$dat.clusters$P
    K <- meta$dat.clusters$K
    
    tbl <- matrix(NA, nrow=0, ncol=N+P)    
    
    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
    
## all (major)
    numEvents <- matrix(NA, nrow=1,ncol=N+P)
    
    gating <- meta$gating
    M  <- rep(0, P)
    inc <- c() 
    
    if( !is.null(gating$childs) ) {
        c <- 0 
        for( g in major ) {
            if( !is.null(gating$childs[[g]]) ) {
                c <- c+1
                
## get mean of major
                for( i in gating$childs[[g]]$clusters ) {
                    j <- which(meta$res.clusters@label==i)
                    M <- M + meta$res.clusters@mu[i,]*
                            sum(meta$dat.clusters$clsEvents[j])
                    
                    inc <- c(inc, j )
                }
                
            }
        }
## mean of major
        M <- M / sum(meta$dat.clusters$clsEvents[inc])
        
        for( p in 1:P ) {
            numEvents[,N+p] <- M[p]
        }
        
## get numEvents of major for all experiment
        k <- 0
        for( i in 1:N ) {
            l <- k+1
            k <- k+K[i]
            ine <- inc[(inc>=l) & (inc<=k)]
            numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine])
        }
        
        for( g in major ) 
        tbl <- rbind(tbl, .meta.numGate(meta, gating$childs[[g]], NULL, 
                                    out.linage=c>1, out.all=out.all))
        
    }
    
    
    parNames <- paste(sep="", "sum(", paste(sep=",", major, collapse=","),")")
    rownames(numEvents) <- parNames
    
    if( out.all ) {
        tbl <- rbind(numEvents, tbl)
    }
    
    for( i in 1:nrow(tbl) ) {
        for( j in 1:N ) {
            tbl[i,j] <- tbl[i,j] / numEvents[1,j]
        }
    }
    
    colnames(tbl) <- c(meta$dat.clusters$expNames, parDesc)
    tbl
}
## meta.relEvents2


#meta.relEvents3 <- function(meta, major=1:5, out.all=T)
#{
#    N <- meta$dat.clusters$N
#    P <- meta$dat.clusters$P
#    K <- meta$dat.clusters$K
#    
#    tbl <- matrix(NA, nrow=0, ncol=N+P)    
#    
#    parDesc <- gsub("\n", ":", meta$res.clusters@parameters)
#    
## all (major)
#    numEvents <- matrix(NA, nrow=(length(major)+1),ncol=N+P)
#    parNames <- c("total", paste(sep="", "P", major))
#    rownames(numEvents) <- parNames
#    
#    numEvents[1,1:N] <- meta$dat.clusters$expEvents
#    
#    gating <- meta$gating
#    M  <- rep(0, P)
#    
#    if( !is.null(gating$childs) ) {
#        majorEvents <- 1
#        for( g in major ) {
#            majorEvents <- majorEvents+1
#            inc <- c()
#            c <- 0
#            if( !is.null(gating$childs[[g]]) ) {
#                c <- c+1
#                
## get mean of major
#                for( i in gating$childs[[g]]$clusters ) {
#                    j <- which(meta$res.clusters@label==i)
#                    M <- M + meta$res.clusters@mu[i,]*
#                             sum(meta$dat.clusters$clsEvents[j])
#                    
#                    inc <- c(inc, j )
#                }
#                
#            }
#            
## mean of major
#            M <- M / sum(meta$dat.clusters$clsEvents[inc])
#            
#            for( p in 1:P ) {
#                numEvents[majorEvents,N+p] <- M[p]
#            }
#            
## get numEvents of major for all experiment
#            k <- 0
#            for( i in 1:N ) {
#                l <- k+1
#                k <- k+K[i]
#                ine <- inc[(inc>=l) & (inc<=k)]
#                numEvents[majorEvents,i] <- 
#                   sum(meta$dat.clusters$clsEvents[ine])
#            }
#            
#            gbl <- .meta.numGate(meta, gating$childs[[g]], NULL, 
#                                 out.linage=c>1, out.all=out.all)
#            for( i in 1:nrow(gbl) ) {
#                for( j in 1:N ) {
#                    gbl[i,j] <- gbl[i,j] / numEvents[majorEvents,j]
#                }
#            }
#            
#            if( out.all ) {
#                for( j in 1:N ) {
#                    numEvents[majorEvents,j] <- numEvents[majorEvents,j] / 
#                                                 numEvents[1,j]
#                }
#                
#                tbl <- rbind(tbl, numEvents[majorEvents,])
#            }
#            
#            rowNames <- rownames(gbl)
#            rownames(gbl) <- paste(sep="/",rowNames, 
#                               rownames(numEvents)[majorEvents])
#            tbl <- rbind(tbl,gbl)
#        }    
#    }
#    
#    
#    colnames(tbl) <- c(meta$dat.clusters$expNames, parDesc)
#    tbl
#}
## meta.relEvents3
### meta.relEvents


###
## meta.parMFI: 
###
.meta.parGate <- 
function(meta, gating, par, name, out.linage=TRUE, out.all=TRUE)
{
    
#    cat("num events\n")
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
        
##        cat(name, out.linage, length(gating$clusters), "\n")
        if( out.linage && length(gating$clusters)>0 ) {
            M  <- rep(0, P)
            inc <- c() 
            for( i in gating$clusters ) {
                j <- which(meta$res.clusters@label==i)
                M <- M + meta$res.clusters@mu[i,]*
                        sum(meta$dat.clusters$clsEvents[j])
                
                inc <- c(inc, j )
            }
            M <- M / sum(meta$dat.clusters$clsEvents[inc])
            
            numEvents <- matrix(NA, nrow=1, ncol=N+P)
            for( p in 1:P ) {
                numEvents[,N+p] <- M[p]
            }
            
## get numEvents for all experiment
            k <- 0
            for( i in 1:N ) {
                l <- k+1
                k <- k+K[i]
                ine <- inc[(inc>=l) & (inc<=k)]
## weighte MFI center
                numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine]*
                                    meta$dat.clusters$M[ine,par])/
                                    sum(meta$dat.clusters$clsEvents[ine])
            }
            
            if( out.all || length(gating$clusters) == 1 ) {
                tbl <- rbind(tbl, numEvents)
                if( length(gating$clusters) == 1 )
                name <- paste(sep=".", name, gating$clusters, 
                            palette()[(gating$clusters%%8)+1])
                rownames(tbl) <- name
            }
        }
        
        if( !is.null(gating$childs) ) {
            c <- 0 
            for( i in 1:length(gating$childs) ) 
            if( !is.null(gating$childs[[i]]) )
            c <- c+1
            
            for( i in 1:length(gating$childs) ) 
            tbl <- rbind(tbl, .meta.parGate(meta, gating$childs[[i]], 
                                        par, name, 
                                        out.linage=c>1, out.all=out.all))
            
        }
        else
        if( length(gating$cluster) > 1 ) {
            rn <- rownames(tbl)
            
            for( cls in gating$cluster ) {
                inc <- which(meta$res.clusters@label==cls)
                numEvents <- matrix(NA, nrow=1, ncol=N+P)                
                k <- 0
                for( i in 1:N ) {
                    l <- k+1
                    k <- k+K[i]
                    ine <- inc[(inc>=l) & (inc<=k)]
## weightes MFI center
                    numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine]*
                                        meta$dat.clusters$M[ine,par])/
                                        sum(meta$dat.clusters$clsEvents[ine])
                }
                for( p in 1:P ) {
                    numEvents[,N+p] <- meta$res.clusters@mu[cls,p]
                }
                
                
                rn <- c(rn, paste(sep=".", name, cls, palette()[cls%%8+1]))
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
    for( i in 1:N ) {
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
    
    tbl <- rbind(tbl, .meta.parGate(meta, meta$gating, par, NULL, 
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
                            palette()[(gating$clusters%%8)+1])
                
                for( p in 1:P ) {
                    numEvents[,N+p] <- meta$res.clusters@mu[gating$clusters,p]
                }
            }
            else {
## weighted sum
                for( p in 1:P ) {
                    numEvents[,N+p] <- M[p]
                }
                
            }
            
## get frequencies for all experiment
            k <- 0
            for( i in 1:N ) {
                l <- k+1
                k <- k+K[i]
                ine <- inc[(inc>=l) & (inc<=k)]
                for( j in 1:nrow(totEvents) ) {
                    numEvents[j,i] <- sum(meta$dat.clusters$clsEvents[ine]) / 
                                        totEvents[j,i]
                }
            }
            
            tbl <- rbind(tbl, numEvents)
            rownames(tbl) <- paste(sep="", name, "./.", rownames(totEvents))
## num events             
            numEvents <- matrix(NA, nrow=1, ncol=N)
            k <- 0
            for( i in 1:N ) {
                l <- k+1
                k <- k+K[i]
                ine <- inc[(inc>=l) & (inc<=k)]
                numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine])
            }
            
            rownames(numEvents) <- name
            parEvents <- rbind(totEvents, numEvents)
            
        }
        
        if( !is.null(gating$childs) ) {
            
            c <- 0 
            for( i in 1:length(gating$childs) ) 
            if( !is.null(gating$childs[[i]]) )
            c <- c+1
            
            for( i in 1:length(gating$childs) ) 
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
                for( i in 1:N ) {
                    l <- k+1
                    k <- k+K[i]
                    ine <- inc[(inc>=l) & (inc<=k)]
                    numEvents[1,i] <- sum(meta$dat.clusters$clsEvents[ine])
                }
                
                for( j in 1:nrow(parEvents) ) {                
                    freqs[j,1:N] <- numEvents[1,] / parEvents[j,]
                    rn <- c(rn, paste(sep=".", name, cls, palette()[cls%%8+1], 
                                    "/", rownames(parEvents)[j]))
                }
                
                for( p in 1:P ) {
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

################################################################################
