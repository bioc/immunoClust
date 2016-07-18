###
## cell event clustering routines
####

###
##  cell.EM
##  fits model to sample data, initial estimation given by K, w, m, s
###
cell.EM <- function(
data, parameters=NULL, expName="immunoClust Experiment", history=NULL, 
state=NULL, K, w, m, s, B=50, tol=1e-5, bias=0.5, modelName="mvt"
) {
    
    y <- .exprs(data, parameters)
    N <- nrow(y)
    P <- ncol(y)
    
    if (nrow(s)) {
        S <- rep(0,length(s))
        for(k in 1:K){
            S[(1+(k-1)*P*P):(k*P*P)] = c(s[k,,])
        }
    }
    else {
        S <- s
    }
    if (nrow(m)) {
        M <- c(t(m))
    }
    else{
        M <- m
    }
    
    obj <- .Call(paste(sep="","immunoC_", modelName, "EMt"), 
                N=as.integer(N), P=as.integer(P), L=as.integer(K), 
                as.double(t(y)), double(0), 
                as.double(w), as.double(M), as.double(S),
                as.integer(B), as.double(tol), as.double(bias) )

    .immunoClust2(obj, K, P, N, state=state,
                    expName=expName, parameters=parameters)
}
### cell.EM

###
##  cell.Estimation
##  classify sample data according to model data given by K, w, m, s
###
cell.Estimation <- function(
data, parameters=NULL, expName="immunoClust Experiment", 
history=NULL, state=NULL, K, w, m, s, modelName="mvt"
) {
    
    y <- .exprs(data, parameters)
    
    y <- as.matrix(y)
    N <- nrow(y)
    P <- ncol(y)
    
    if (nrow(s)) {
        S <- rep(0,length(s))
        for(k in 1:K){
            S[(1+(k-1)*P*P):(k*P*P)] = c(s[k,,])
        }
    }
    else {
        S <- s
    }
    if (nrow(m)) {
        M <- c(t(m))     
    }
    else{
        M <- m
    }
    
    obj <- .Call(paste(sep="", "immunoC_", modelName, "E"), 
                as.integer(N), as.integer(P), as.integer(K), 
                as.double(t(y)), double(0), 
                as.double(w), as.double(M), as.double(S) )

    .immunoClust2(obj, K, P, N, expName=expName, parameters=parameters);
    
}
### cell.Estimation


### cell.ME
###
##  fit model to the sample data, initial event assignment given by label
###
cell.ME<-function(
data, parameters=NULL, expName="immunoClust Experiment", 
history=NULL, state=NULL, label, B=50, tol=1e-5, modelName="mvt"
) {
    
    y <- .exprs(data, parameters);
    
    inc <- !is.na(label)
    
    y <- as.matrix(y[inc,])
    N <- nrow(y)
    P <- ncol(y)
    label <- label[inc]
    K <- max(label)
    
    obj <- .Call(paste(sep="", "immunoC_", modelName, "ME"), 
            as.integer(N), as.integer(P), as.integer(K), 
            as.double(t(y)), NULL,  as.integer(label),
            as.integer(B), as.double(tol) )

    .immunoClust2(obj, K, P, N, expName=expName, parameters=parameters, inc=inc)
}
### cell.ME


### cell.FitModel
###
##  fit pre-assumed model x to sample data
###
cell.FitModel <- 
function(x, data, B=50, tol=1e-5, bias=0.5, modelName="mvt" ) 
{
    s <- strptime(date(), "%a %b %d %H:%M:%S %Y")
    
    res <- cell.EM(data, parameters=x@parameters, 
                history=attr(x, "history"), state=attr(x,"state"),
                K=x@K, w=x@w, m=x@mu, s=x@sigma, 
                B=B, tol=tol, bias=bias, modelName=modelName)
    
    attr(res, "trans.a") <- attr(x,"trans.a")
    attr(res, "trans.b") <- attr(x,"trans.b")
    attr(res, "trans.decade") <- attr(x,"trans.decade")
    attr(res, "trans.scale") <- attr(x,"trans.scale")
    
    e <- strptime(date(), "%a %b %d %H:%M:%S %Y")
    message("EM takes ", format(difftime(e,s,units="min"), digits=2), 
            " minutes\n")
    
    res
}
### cell.FitModel

### cell.Classify
###
##  assign sample data to model
###
cell.Classify <- function(x, data, modelName="mvt" ) {
    
    res <- cell.Estimation(data, parameters=x@parameters, 
                    history=x@history, state=attr(x,"state"),
                    K=x@K, w=x@w, m=x@mu, s=x@sigma, modelName=modelName)
    
    attr(res, "trans.a") <- attr(x,"trans.a")
    attr(res, "trans.b") <- attr(x,"trans.b")
    attr(res, "trans.decade") <- attr(x,"trans.decade")
    attr(res, "trans.scale") <- attr(x,"trans.scale")
    
    res
}
### cell.Classify


### cell.FitCluster
###
##  carry cluster-asignment to sample data 
##  (especially uncompensated model/assignment to compensated sample data)
###
#cell.FitCluster <- function(x, data, B=50, tol=1e-5, modelName="mvt" ) {
# res <- cell.ME(data, parameters=x@parameters, 
#                history=attr(x, "history"), state=attr(x, "state"),
#                label=x@label, B=B, tol=tol, modelName=modelName)
# 
# attr(res, "trans.a") <- attr(x, "trans.a")
# attr(res, "trans.b") <- attr(x, "trans.b")
# attr(res, "trans.decade") <- attr(x, "trans.decade")
# attr(res, "trans.scale") <- attr(x,"trans.scale")
#
# res 
#}
### cell.FitCluster

### cell.ClustData
##  cluster sample data in K clusters
###
cell.ClustData<-function(
data, K, parameters=NULL, expName="immunoClust Experiment", 
sample.seed=1, sample.number=1500, sample.standardize=TRUE,
B=50, tol=1e-5, modelName="mvt"
) {
    
    y <- .exprs(data, parameters)
    
    N <- nrow(y)
    P <- ncol(y)
    
    
# to perform the cluster analysis via EM for each specific number of clusters
    if (K==1) { 
        label <- rep(1, N)
    }
    else {
        if ( (P==1) ) {
            q <- quantile(y, seq(from=0, to=1, by=1/K))
            label <- rep(0, N)
            q[1] <- q[1]-1
            for (k in 1:K) label[y>q[k] & y<=q[k+1]] <- k
        }
        else {
            if (N > sample.number) {
                set.seed(sample.seed)
                ySubset <- sample(1:N, sample.number)
            }
            else {
                ySubset <- 1:N
            }
            
## 2013.01.30: does not matter
            if( sample.standardize ) {
                x <- y[ySubset,]
                for( p in 1:P )
                x[,p] <- (x[,p]-mean(x[,p]))/sd(x[,p])
            }
            
            hcPairs <- cell.hclust(y[ySubset,])
            
            label <- rep(0, N)
            label[ySubset] <- .clust.hclass(hcPairs, K)
        }
    }
    
# EMs
    obj <- .Call(paste(sep="", "immunoC_", modelName, "ME"), 
                as.integer(N), as.integer(P), as.integer(K), 
                as.double(t(y)), NULL, as.integer(label), 
                as.integer(B), as.double(tol) )
    
    .immunoClust2(obj, K, P, N, expName=expName, parameters=colnames(y))
}
### cell.ClustData


###
### cell.SubClustering
###
##  try sub-clustering for each cluster and select most increasing 
##  sub-clustering for the whole model
###

###
##  replace cluster k in x by clusters in y
###
.mergeModel <- function( x, y, k)
{
    K <- x@K
    P <- dim(x@mu)[2]
    L <- y@K
    
    W <- rep(0, (K+L-1))
    
    history <- rep("", (K+L-1))
    state <- rep(0,(K+L-1))
    
    M <- rep(0, (K+L-1)*P)
    dim(M) <- c(K+L-1,P)
    S <- rep(0, (K+L-1)*P*P)
    dim(S) <- c(K+L-1, P, P)
    
    if( k>1 ) {
        for( i in 1:(k-1) ){ 
            W[i] = x@w[i]
            M[i,] = x@mu[i,]
            S[i,,] = x@sigma[i,,]
            if( !is.null(x@history) ) {
                history[i] <- x@history[i]
            }
            if( !is.null(x@state) ) {
                state[i] <- x@state[i]
            }
        }
    }
    if( k<K ) {
        for( i in (k+1):K ) {
            W[L+i-1] = x@w[i]
            M[L+i-1,] = x@mu[i,]
            S[L+i-1,,] = x@sigma[i,,]
            if( !is.null(x@history) ) {
                history[L+i-1] <- x@history[i]
            }
            if( !is.null(x@state) ) {
                state[L+i-1] <- x@state[i]
            }
        }
    }
    
    for( i in 1:L ){
        W[k-1+i] = x@w[k]*y@w[i]
        S[k-1+i,,] = y@sigma[i,,]
        M[k-1+i,] = y@mu[i,]
        if( !is.null(x@history) ) {
            history[k-1+i] <- x@history[k]
        }
    }
    res <- new("immunoClust", expName="Model Refinement", 
            parameters=x@parameters, 
            K=K+L-1, N=x@N, P=x@P ,w=W, mu=M, sigma=S, 
            history=history, state=state)
    res
}


cell.SubClustering <- function(
x, dat, B=50, tol=1e-5, thres=0.1, bias=0.5,
sample.weights=1, sample.EM="MEt", sample.number=1500, 
sample.standardize=TRUE, extract.thres=0.8, modelName="mvt"
) {
    
    s <- strptime(date(), "%a %b %d %H:%M:%S %Y")
    
    y <- .exprs(dat, x@parameters)
    
    inc <- !is.na(x@label)  
    y <- as.matrix(y[inc,])
    z <- as.matrix(x@z[inc,]) 
    inc <- inc[inc] 
    
## 2014.05.07: ungluecklicker quick fix
    if( sum(is.infinite(z) | is.na(z) | is.nan(z)) > 0 ) {
        warning("Fehler: Z has infinite values", 
                sum(is.infinite(z) | is.na(z) | is.nan(z) ), "\n")
        z[is.infinite(z) | is.na(z) | is.nan(z)] <- 0
    }
    
    N <- nrow(y)
    P <- ncol(y)
    K <- x@K
    
    state <- x@state
    if( is.null(state) || length(state) != x@K ) {
        state <- rep(numeric(0), x@K)
    }
    
## test sub-models with 1 to 8 sub-clusters
    J <- 8
    cutoff <- 0
    
    icl_thres <- (P*(P+1)/2 + P)*log(N)*0.5*thres
    icl_OK <- (P*(P+1)/2 + P)*log(N)*0.5
    
    model <- new("immunoClust", expName="Model Refinement", 
                parameters=x@parameters,
                N=N, P=P, K=x@K, w=x@w, mu=x@mu, sigma=x@sigma, 
                history=x@history, state=x@state)
    
    res_l <- vector("list", K)
    icl_l <- rep(0, K)
    tst_l <- rep(1, K)
    
    state = model@state
    for( k in 1:K ) {
        
        message("Test cluster ", k, " for sub-clustering")
## get cluster data
        cinc <- .clusterData(y,z,inc, k, extract.thres)
        t <- NULL 
## 2014.06.16: use weights T?
        w <- as.double(sample.weights)
        if( is.double(w) & w > 0 ) {  
            t <- z[cinc,k]^w
        }   
        
        ks <- strptime(date(), "%a %b %d %H:%M:%S %Y")
        
        res <- cell.TestSubCluster( x, as.matrix(y[cinc,]), t, k, J=J, 
                            B=B, tol=tol, bias=bias,
                            sample.EM=sample.EM, sample.number=sample.number, 
                            sample.standardize=sample.standardize, 
                            modelName=modelName) 
        
        ke <- strptime(date(), "%a %b %d %H:%M:%S %Y")
        message("EM takes ", format(difftime(ke,ks,units="min"), digits=2), 
                " minutes")

#cat("cluster", k, "(", length(cinc), ")", "state", model@state[k])
        
        res_l[[k]] <- res
        if( !is.null(res) && length(res) > 1 ) {
        
            icl <- rep(0, length(res)-1)

            for( l in 2:length(res) )
            icl[l-1] <- res[[l]]@ICL/res[[l]]@K
            
            icl_l[k] <- max(icl)
            l <- 1+which.max(icl)
            tst_l[k] <- l

            
            model@state[k] <- icl_l[k]
#           if( icl_l[k] < -icl_OK && res[[l]]@K==1 ) {
#               model@state[k] <- 5
#           }
#           else
#           if( icl_l[k]*res[[l]]@K < -icl_OK ) {
#               model@state[k] <- 5
#           }
#           else
#           if( icl_l[k] < -icl_OK ) {
#               model@state[k] <- 4
#           }
#           else
#           if( res[[l]]@K == 1 ) {
#               model@state[k] <- 3
#           }
#           else
#           if( icl_l[k] < icl_thres ) {
#               model@state[k] <- 1
#           }
#           else {
#               model@state[k] <- 0
#           }
            
#if( res[[2]]@ICL < -2*bias ) {
#cat("cluster", k, "state", model@state[k], "=> 1\n")
#               model@state[k] <- 1
#            }
        }
        else {
## state already > 0
#     model@state[k] <- 2
            icl_l[k] <- 0
            tst_l[k] <- 1
        }

#cat("=>", model@state[k], "\n")
        
    } ## for cluster k
    
    
    off <- 0
    ins <- vector("list",K)
    sK <- 0
    xK <- max(8,2*x@K)
    
    while( xK > sK ) {
        k <- which.max(icl_l)
        
        res <- res_l[[k]]
        icl <- icl_l[k]
        l <- tst_l[k]
        
        if( is.null(res) ) {
            break
        } 
        
        if( res[[l]]@K > 1 ) {
            message("cluster ", k, " has ", res[[l]]@K, " sub-cluster at ", l, 
                    ", ICL=", format(icl, digits=2))
        }
        
        icl_l[k] <- cutoff
        
        res <- res[[l]]
        
        if( icl <= cutoff )
        break
        
        if( (res@K>1) & (icl>icl_thres) ) {
            ins[[k]] <- new("immunoClust", expName="Cluster Refinement", 
                        parameters=res@parameters,
                        K=res@K, w=res@w, mu=res@mu, sigma=res@sigma, 
                        state=rep(0, res@K) )
        }   
        
        sK <- 0
        for(i in 1:K ) if( !is.null(ins[[i]]) ) 
        sK <- sK + (ins[[i]])@K
    } 
    
    
    for( k in 1:K) if( !is.null(ins[[k]]) ) {
#        message("split cluster ", k, " into ", (ins[[k]])@K, " sub-cluster")
#cat("split cluster", k, "state:", state[k], "=>", model@state[k+off], "\n")
        model <- .mergeModel(model, ins[[k]], k+off)
        off <- off + (ins[[k]]@K) - 1
    }  
    
    attr(model, "trans.a") <- attr(x,"trans.a")
    attr(model, "trans.b") <- attr(x,"trans.b")
    attr(model, "trans.decade") <- attr(x,"trans.decade")
    attr(model, "trans.scale") <- attr(x,"trans.scale")
    
    e <- strptime(date(), "%a %b %d %H:%M:%S %Y")
    message("Model Refinement takes ", 
            format(difftime(e,s,units="min"), digits=2), " minutes\n")
    
    model
}
### cell.SubClustering

### cell.TestSubCluster
###
##  called by cell.SubClustering, 
##  calculates clustering for 1 to J cluster on sample flowFrame Y
###
.icl_lambda <- function(N,P,L) {
    (L-1)*(P*(P+1)/2 + P)*log(N)*0.5
}
.icl_delta <- function(N,P,K, L) {
## delta for L sub cluster if one cluster
    res <- (L-1)*(P*(P+1)/2 + P)*log(N)*0.5
    res <- res - (lgamma((K+L-1)/2) - lgamma(K/2))
    res <- res + (L-1) * lgamma(1/2)
    res <- res + lgamma(N+(K+L-1)/2) - lgamma(N+K/2)
    res
}
.icl_delta2 <- function(N,K,L) {
    res <- 0.0
    res <- res - (lgamma((K+L-1)/2) - lgamma(K/2))
    res <- res + (L-1) * lgamma(1/2)
    res <- res + lgamma(N+(K+L-1)/2) - lgamma(N+K/2)
    res
}

.clusterInclude <- function(x, y, inc, cluster, thres=0.99)
{
    P <- ncol(y)
    N <- nrow(y)
    K = x@K
    
    if (nrow(x@sigma)) {
        S <- rep(0,length(x@sigma))
        for(k in 1:K){
            S[(1+(k-1)*P*P):(k*P*P)] = c(x@sigma[k,,])
        }
    }
    else {
        S <- x@sigma
    }
    if (nrow(x@mu)) {
        M <- c(t(x@mu))    
    }
    else{
        M <- x@mu
    }
    
    ret <- .Call("immunoC_clusterInclude",
                as.integer(N), as.integer(P), as.integer(K),
                as.double(t(y)), as.double(x@w), as.double(M), as.double(S),
                as.integer(cluster), as.integer(inc), as.double(thres))
    
    which(as.logical(ret))
}

.clusterData <- function(y, z, inc, cluster, thres=0.8)
{
    P <- ncol(y)
    N <- nrow(y)
    K <- ncol(z)
    ret <- .Call("immunoC_clusterData",
                as.integer(N), as.integer(P), as.integer(K),
                as.double(NULL), as.double(t(z)),
                as.integer(cluster), as.integer(inc), as.double(thres))
    
    which(as.logical(ret))
}

cell.TestSubCluster<-function(
x, y, t, cluster, J=8, B=500, tol=1e-5, bias=0.5,
sample.EM="MEt", sample.df=5, sample.number=1500, sample.standardize=TRUE, 
modelName="mvt"
) {
## total model
    N <- nrow(y)
    P <- ncol(y)
    K <- x@K
    
    sumT <- N
    if( !is.null(t) ) 
    sumT <- sum(t)
    
    tY <- t(y)
    
    prob <- NULL  
    
    if( J > N ) {
        return(NULL)
    }
    
    result <- vector("list", J)
    
    label <- rep(1, N)
    
    obj <- .Call(paste(sep="", "immunoC_", modelName, "ME"), 
                as.integer(N), as.integer(P), L=as.integer(1), 
                as.double(tY), as.double(t), as.integer(label),
                as.integer(B), as.double(tol))  
    
    
    if( obj$L < 1 ) 
    return(NULL)
    
# output obj$s to sigma
    sigma <- array(0, c(1, P, P))
    s <- matrix(obj$s, 1, P * P, byrow=TRUE)
    sigma[1,,] <- matrix(s[1,], P, P, byrow = TRUE)
    
    
# output BIC & ICL
    BIC <- obj$logLike[1]
    ICL <- 0
    logLike <- obj$logLike[3]
    iclLike <- obj$logLike[2]
    
# outp    
    result[[1]] <- new("immunoClust", parameters=x@parameters, 
                    K=1, N=N, P=P, w=obj$w, 
                    mu=matrix(obj$m, 1, P, byrow=TRUE), sigma=sigma,
                    logLike=obj$logLike, BIC=BIC, ICL=ICL)   
    
    obj <- NULL 
# initialization based on hierarchical clustering
    if (J>1 && P>1 ) {
        prob <- NULL
        maha <- NULL
        
        use_p <- which(diag(x@sigma[cluster,,]) > 1e-8)
        try( maha <- mahalanobis(y[,use_p], x@mu[cluster,use_p], 
                            x@sigma[cluster, use_p, use_p]) )
        
        if( is.null(maha) ) {
#warning(" singularity:", det(x@sigma[cluster,,]), 
#                   diag(x@sigma[cluster,,]), "use", use_p, "\n") 
            warning(" singularity in cluster", cluster, "\n")
        }
        else {
            abv <- qchisq(0.95,P)^2
            maha[maha>abv] <- abv
            
# density based down sampling: should be adapted to model (mvt -> t, mvn -> n)
            if( "mvt" == modelName ) {
                prob <- (1 + maha/sample.df)^(0.5*(sample.df+P))
            }
            else
            if( "mvn" == modelName ) {
                prob <- exp(0.5*maha)
            }
        }
        
# if more than sample.number (1500) observations, only use testSample at random
        if (N > sample.number) {
# 2012.11.07:     
# enhence outliers in sub sample
            if( !is.null(maha) ) {
                ySubset <- sample(1:N, sample.number, prob=prob)
            }
            else {
                ySubset <- sample(1:N, sample.number)
            }
        }
        else {
            ySubset <- 1:N
        }
        
##    hcPairs <- HClust(y[ySubset,], weights=t[ySubset])
## ... or standardize?      
        sub <- y[ySubset,]
        if( sample.standardize ) {
            for( p in 1:P ) {
                if( sd(sub[,p]) > 0 ) {
                    sub[,p] <- (sub[,p]-mean(sub[,p]))/sd(sub[,p])
                }
            }
        }
        
        hcPairs <- cell.hclust(sub, t[ySubset])
        
        attr(hcPairs, "ySubset") <- ySubset
    }
    
##  to perform the cluster analysis via EM for each specific number of clusters
    if( J > 1 ) for (k in 2:J) {
        
        obj <- NULL
        gc(verbose=FALSE, reset=TRUE)
        
        label <- rep(0, N)
        
## TODO: clarify P=1 case   
        if (P==1) {
            q <- quantile(y, seq(from=0, to=1, by=1/k))
            q[1] <- q[1]-1
            for (l in 1:k) label[y>q[l] & y<=q[l+1]] <- l
        }
        else {
            label[ySubset] <- .clust.hclass(hcPairs, k)
        }   
        
# EMs        
        obj <- .Call(paste(sep="", "immunoC_", modelName, sample.EM), 
                    as.integer(N), as.integer(P), as.integer(k), 
                    as.double(tY), as.double(t), as.integer(label), 
                    as.integer(B), as.double(tol), as.double(bias) ) 
        
## 2012.12.12: singularity problems   
        if( obj$L < 1 || obj$logLike[3] == Inf || obj$tolerance > tol) {
            res_t = vector("list", k-1)
            for( l in 1:(k-1) )
            res_t[[l]] <- result[[l]]
            
            result <- res_t
            J = k-1
            break
        }
        
        L <- obj$L
# output obj$s to sigma
        sigma <- array(0, c(L, P, P))
        s <- matrix(obj$s, k, P * P, byrow=TRUE)
        for (l in 1:L)
        sigma[l,,] <- matrix(s[l,], P, P, byrow = TRUE)
        
        mu <- matrix(obj$m, k, P, byrow=TRUE)[1:L,]
        dim(mu) <- c(L,P)
        
# output BIC & ICL
        BIC <- obj$logLike[1]
        
## 2012.11.07: use sumT not total N
## 2012.12.13: use obj$L and not k  
        ICL <- obj$logLike[3] - logLike - .icl_delta(sumT, P, K, L)*bias
## 2016.06.28: skip below, is a bit unpredictable      
#        if( L > result[[k-1]]@K ) {
#            DCL <- obj$logLike[3] - result[[k-1]]@logLike[3] - 
#                    .icl_delta(sumT, P, K, L)*bias
#            if( DCL > 0 && DCL > ICL ) {
#                ICL <- DCL
#            }
#        }
        
# outp    
        result[[k]] <- new("immunoClust", parameters=x@parameters, 
                        K=L, N=N, P=P, w=obj$w[1:L], mu=mu, sigma=sigma,
                        logLike=obj$logLike, BIC=BIC, ICL=ICL)
        obj <- NULL
        
    } ## for k
    
    result
    
}
### cell.TestSubCluster

###
##  reimplementation of hcvvv in mclust with optional weights
###
cell.hclust <- function(data, weights=NULL)
{
    
    if(any(is.na(data)))
    stop("missing values not allowed in data")
    
    data <- as.matrix(data)
    dimdat <- dim(data)
    
    if(is.null(dimdat) || length(dimdat) > 2)
    stop("data should in the form of a matrix")
    
    dimnames(data) <- NULL
    N <- nrow(data)
    P <- ncol(data)
    
    if(N <= P)
    warning("# of observations <= # of parameters")
    
    partition <- 1:N
    attr(partition, "unique") <- N
        
    obj <- .Call("immunoC_mvnHC", as.integer(N), as.integer(P), 
                as.double(t(data)), as.double(weights))
    
    structure(t(cbind(obj$li,obj$lj)), 
            initialPatition=partition, change=obj$crit,
            dimensions=c(N,P), modelName = "mvn", 
            call = match.call())
}
### cell.hclust
