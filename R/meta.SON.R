
meta.SON.preclustering <- function(
exp, par.subset=c(), meta.iter=5, tol=1e-3, meta.bias=1,
meta.alpha=.5, HC.samples= 500, verbose=FALSE
)
{
    dat <- meta.exprs(exp, sub=par.subset)
    
    res <- meta.Clustering(dat$P, dat$N, dat$K,
                        dat$clsEvents, dat$M, dat$S,
                        I.iter=meta.iter, B=10, tol=tol,
                        bias=meta.bias, sub.thres=meta.bias, alpha=meta.alpha,
                        HC.samples=HC.samples,
                        verbose=verbose)

    
    attr(res, "desc") <- dat$desc
    if( is.null(par.subset) )
        par.subset <- seq_len(npar(res))
    attr(res, "trans.a") <- apply( vapply(exp,function(x) x@trans.a[par.subset],
                                rep(0.01,npar(res))), 1, mean)
    attr(res, "trans.b") <- apply( vapply(exp,function(x) x@trans.b[par.subset],
                                rep(0.0,npar(res))), 1, mean)
    attr(res, "limits") <- attr(exp[[1]], "limits")[,par.subset]

    meta <- list("dat.scatter"=NULL, "res.scatter"=NULL,
                    "dat.clusters"=dat, "res.clusters"=res)
    meta$gating <- list("clusters"=seq_len(res@K), "childs"=c(),
                    "desc"="all", "partition"=TRUE)

    
    meta$gating$pscales <- Default_Scales(attr(res, "trans.a"),
                                        attr(res, "limits"))
    
    class(meta) <- "immunoMeta"
    meta
}

##
# meta.SON.clustering
##
#
meta.SON.clustering <- function(
meta,
cycles=6, alpha=0.5, scale.factor=2, scale.steps=0,
meta.iter=1, meta.bias=0.3, meta.thres=meta.bias, meta.tol=1e-5,
SON.cycles=1, SON.rlen=100, SON.deltas=c(1/SON.rlen,1/SON.rlen),
SON.blurring=c(2,0.1),
verbose=FALSE
)
{
    dat <- meta$dat.clusters
    res <- meta$res.clusters
    ndat <- dat
    nres <- res
    
    nM <- dat$M
    sub.iter <- 0
    for( cycle in seq_len(cycles) ) {
        if(verbose) {
            #message("SON.clustering cycle", cycle, "\n")
            message("SON/ormalize", cycle, "\n")
        }
        ## allways original M or normedM?
        s <- strptime(date(), "%a %b %d %H:%M:%S %Y")
        
        attr(nres,"P") <- dat$P
        obj <- .Call("immunoC_SON_normalize",
            nres, as.integer(dat$N), as.integer(dat$K),
            as.double(dat$clsEvents), 
            as.double(t(dat$M)),
            ##as.double(t(nM)),
            as.double(t(dat$S)),
            as.double(alpha), as.double(scale.factor), as.integer(scale.steps),
            as.integer(sub.iter), as.double(meta.tol),
            as.integer(SON.cycles), as.integer(SON.rlen), as.double(SON.deltas),
            as.double(SON.blurring)
            )
        
        
        nM <- matrix(obj, nrow(dat$M), ncol(dat$M), byrow=TRUE)
        colnames(nM) <- colnames(dat$M)
        e <- strptime(date(), "%a %b %d %H:%M:%S %Y")
        
        if(verbose) {
            message( "\ttakes ", format(difftime(e,s,units="min")), "\n")
            message("SON Clustering", cycle, "\n")
        }
        s <- e
        nres <- meta.Clustering(dat$P, dat$N, dat$K, dat$clsEvents, nM, dat$S,
            label=nres@label, I.iter=meta.iter, tol=meta.tol,
            bias=meta.bias, sub.thres=meta.thres, alpha=alpha,
            EM.method=200,
            ## HC.samples?
            verbose=verbose)
        
        e <- strptime(date(), "%a %b %d %H:%M:%S %Y")
        
        if(verbose)
        message("\ttakes ", format(difftime(e,s,units="min")), "\n" )
    }
    ## 2024.04.17: will man das so?
    dat$nrm.M <- nM
    attr(nres,"trans.a") <- attr(res, "trans.a")
    attr(nres,"trans.b") <- attr(res, "trans.b")
    attr(nres,"limits") <- attr(res, "limits")
    ret <- immunoMeta(nres, dat)
    ##  2024.04.17: oder eigentlich so?
    ##ret$dat.clusters$nrm.M <- nM
    parameters(ret) <- parameters(meta)
    prop(ret,"pscales",c()) <- prop(meta,"pscales",c())
    attr(res,"SON.call") <- match.call()
    
    ret
    
}

meta.SON.combineClustering <-
function(meta, res, par=seq_len(npar(meta)),
map.cluster=seq_len(ncls(meta)), use.cluster=seq_len(ncls(res)),
meta.alpha=0.5, meta.bias=0.1, meta.iter=100, meta.tol=1e-5,
SON.method=1, SON.cycles=4, SON.rlen=10,
SON.deltas=c(1/SON.rlen,1/SON.rlen), SON.blurring=c(2,1),
traceG=c(), traceK=c())
{
    ## co-clustering with normalization
    model_res <- subset(meta$res.clusters, par=par)
    model_clsEvents <- meta$dat.clusters$clsEvents
    
    sample_res <- subset(res, par=par)
    sample_clsEvents <- events(res)
    # equalize the numberof total events
    sample_clsEvents <- sample_clsEvents *
    sum(model_clsEvents)/sum(sample_clsEvents)
    # build context 20:1 model/test events
    ## models should be normalized to N=1: weight src
    ## dst by 10:1 maybe a parameter?
    sample_clsEvents <- sample_clsEvents / 20
    
    attr(model_res, "evts") <- as.double(model_clsEvents)
    attr(sample_res,"evts") <- as.double(sample_clsEvents)
    
    
    map_cluster <- rep(0, ncls(model_res))
    use_cluster <- rep(0, ncls(sample_res))
    map_cluster[map.cluster] <- 1
    use_cluster[use.cluster] <- 1
    
    obj <- .Call("immunoC_SON_combineClustering",
    model_res, sample_res, as.integer(map_cluster), as.integer(use_cluster),
    as.double(meta.alpha),  as.double(meta.bias),
    as.integer(meta.iter), as.double(meta.tol),
    as.integer(SON.method), as.integer(SON.cycles), as.integer(SON.rlen),
    as.double(SON.deltas), as.double(SON.blurring),
    as.integer(c(traceG,0)), as.integer(c(traceK,0))
    )
    
    co_res <- .immunoClust2(obj, ncls(model_res), npar(model_res),
    ncls(model_res)+ncls(sample_res),
    parameters=parameters(model_res))
    
    normedM <- matrix(obj$normedM, nobs(co_res), npar(co_res), byrow=TRUE)
    
    attr(co_res, "nrm.M") <- normedM
    
    co_res
}

meta.SON.normalize <- function(
meta,
alpha=0.5, scale.factor=2, scale.steps=0,
SON.cycles=1, SON.rlen=100, SON.deltas=c(1/SON.rlen,1/SON.rlen),
SON.blurring=c(2,0.1),
verbose=FALSE
)
{
    dat <- meta$dat.clusters
    res <- meta$res.clusters
    ndat <- dat
    nres <- res
    
    nM <- dat$M
    sub.iter <- 0
    
    #for( cycle in seq_len(cycles) ) {
       
        ## allways original M or normedM?
        s <- strptime(date(), "%a %b %d %H:%M:%S %Y")
        
        attr(nres,"P") <- dat$P
        obj <- .Call("immunoC_SON_normalize",
            nres, as.integer(dat$N), as.integer(dat$K),
            as.double(dat$clsEvents),
            as.double(t(dat$M)),
            ##as.double(t(nM)),
            as.double(t(dat$S)),
            as.double(alpha), as.double(scale.factor), as.integer(scale.steps),
            as.integer(0), as.double(0.0), ## unused
            as.integer(SON.cycles), as.integer(SON.rlen), as.double(SON.deltas),
            as.double(SON.blurring)
            )
        
        
        nM <- matrix(obj, nrow(dat$M), ncol(dat$M), byrow=TRUE)
        colnames(nM) <- colnames(dat$M)
        e <- strptime(date(), "%a %b %d %H:%M:%S %Y")
        
    #}
    
    nM
    
}
