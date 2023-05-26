

##
# meta.SON.clustering
##
#
meta.SON.clustering <- function(
meta,
cycles=6, alpha=0.5, scale.factor=2, scale.steps=0,
meta.iter=2, meta.bias=0.3, meta.tol=1e-5,
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
        if(verbose)
        message("SON.clustering cycle", cycle, "\n")
        ## allways original M or normedM
        
        attr(nres,"P") <- dat$P
        obj <- .Call("immunoC_SON_normalize",
            nres, as.integer(dat$N), as.integer(dat$K),
            as.double(dat$clsEvents), as.double(t(dat$M)), as.double(t(dat$S)),
            as.double(alpha), as.double(scale.factor), as.integer(scale.steps),
            as.integer(sub.iter), as.double(meta.tol),
            as.integer(SON.cycles), as.integer(SON.rlen), as.double(SON.deltas),
            as.double(SON.blurring)
            )
        
        
        nM <- matrix(obj, nrow(dat$M), ncol(dat$M), byrow=TRUE)
        colnames(nM) <- colnames(dat$M)
        
        nres <- meta.Clustering(dat$P, dat$N, dat$K, dat$clsEvents, nM, dat$S,
            label=nres@label, I.iter=meta.iter, tol=meta.tol,
            bias=meta.bias, sub.thres=1, alpha=alpha, norm.method=0,
            verbose=verbose)
        
        if(verbose)
        message("\t=>", ncls(nres), "clusters\n" )
    }
    dat$M <- nM
    attr(nres,"trans.a") <- attr(res, "trans.a")
    attr(nres,"trans.b") <- attr(res, "trans.b")
    attr(nres,"limits") <- attr(res, "limits")
    ret <- immunoMeta(nres, dat)
    #ret$dat.clusters$nrm.M <- nrm_M
    parameters(ret) <- parameters(meta)
    prop(ret,"pscales",c()) <- prop(meta,"pscales",c())
    
    ret
    
}

meta.SON.combineClustering <-
function(meta, res, par=seq_len(npar(meta)),
map.cluster=seq_len(ncls(meta)), use.cluster=seq_len(ncls(res)),
meta.alpha=0.5, meta.cycles=1, meta.bias=0.1,
meta.iter=100, meta.tol=1e-5,
SON.cycles=4, SON.rlen=10,
SON.deltas=c(1/SON.rlen,1/SON.rlen), SON.blurring=c(2,1), SON.norm=1,
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
    as.double(meta.alpha), as.integer(meta.cycles), as.double(meta.bias),
    as.integer(meta.iter), as.double(meta.tol),
    as.integer(SON.cycles), as.integer(SON.rlen), as.double(SON.deltas),
    as.double(SON.blurring), as.integer(SON.norm),
    as.integer(c(traceG,0)), as.integer(c(traceK,0))
    )
    
    co_res <- .immunoClust2(obj, ncls(model_res), npar(model_res),
    ncls(model_res)+ncls(sample_res),
    parameters=parameters(model_res))
    
    normedM <- matrix(obj$normedM, nobs(co_res), npar(co_res), byrow=TRUE)
    
    attr(co_res, "nrm.M") <- normedM
    
    co_res
}

