setClass("immunoClust",
representation(expName="character", fcsName="character", parameters="character",
removed.below="numeric", removed.above="numeric",
trans.a="numeric", trans.b="numeric", 
trans.decade="numeric", trans.scale="numeric",
K="numeric", N="numeric", P="numeric", w="vector", mu="matrix", sigma="array",  
z="matrix", label="vector", 
logLike="numeric", BIC="numeric", ICL="numeric", 
history="character", state="vector"
),
prototype(expName=character(0), fcsName=character(0), parameters=character(0),
removed.below=numeric(0), removed.above=numeric(0), 
trans.a=numeric(0), trans.b=numeric(0), trans.decade=-1, trans.scale=1,  
K=0, N=0, P=0, w=rep(numeric(0),0), mu=matrix(numeric(0), nrow=0, ncol=0),
sigma=array(numeric(0), c(0,0,0)), 
z=matrix(numeric(0), nrow=0, ncol=0), label=rep(numeric(0),0),
logLike=numeric(0), BIC=numeric(0), ICL=numeric(0), 
history=character(0), state=rep(numeric(0),0)
)
)

## show method

setMethod("show", "immunoClust",
function(object)
{
    cat("Object of class 'immunoClust'","\n")
    cat("This object has the following slots: \n")
    cat("expName, fcsName ,parameters,",
        "removed.below, removed.above,",
        "trans.a, trans.b, trans.decade, trans.scale",
        "K, N, P, w, mu, sigma, z, label,",
        "logLike, BIC, ICL\n")
})


## summary method
setMethod("summary", "immunoClust",
function(object)
{
    param <- attributes(object)$param
    P <- length(param)
    N <- length(object@label)
    cat("** Experiment Information ** \n")
    cat("Experiment name:", attr(object, "expName"),"\n")
    cat("Data Filename:  ", attr(object, "fcsName"),"\n")
    cat("Parameters:  ", param,"\n")  
    cat("Description: ", attr(object, "desc"), "\n")
    
    cat("\n** Data Information ** \n")
    cat("Number of observations:",N,"\n")
    cat("Number of parameters:  ",P,"\n")
    if( length(object@removed.above) > 0 && length(object@removed.below) > 0 ){
    cat("Removed from above:    ", object@removed.above, 
        " (", round(object@removed.above/N*100, 2), "%)\n", sep="")

    cat("Removed from below:    ", object@removed.below,
        " (", round(object@removed.below/N*100, 2), "%)\n", sep="")
    }
    else {
    removed <- sum(is.na(object@label))
    cat("Removed observations:   ", removed,
        " (", round(removed/N*100, 2), "%)\n", sep="")
    }
    
    cat("\n** Transformation Information ** \n")
    cat("htrans-A:  ",formatC(attr(object,"trans.a"),format="f",digits=6),"\n")
    cat("htrans-B:  ",formatC(attr(object,"trans.b"),format="f",digits=6),"\n")
    cat("htrans-decade:  ", attr(object, "trans.decade"),"\n")
    
    cat("\n** Clustering Summary ** \n")
    cat("ICL bias:", formatC(attr(object,"bias"),format="f",digits=2),"\n")
    cat("Number of clusters:",object@K,"\n")
    
    cat("Cluster     Proportion  Observations\n")
    for( k in seq_len(object@K) ) {
        cat( sprintf("%8d    %10.6f  %12d\n", k, object@w[k], 
            sum(!is.na(object@label) & object@label==k) ) )
    }   
    cat("\n")
    k <- which.min(object@w)
    cat( sprintf("%8s    %10.6f  %12d\n", "Min.", object@w[k], 
        sum(!is.na(object@label) & object@label==k) ) )
    k <- which.max(object@w)
    cat( sprintf("%8s    %10.6f  %12d\n", "Max.", object@w[k], 
        sum(!is.na(object@label) & object@label==k) ) )
    
    cat("\n** Information Criteria ** \n")
    cat("Log likelihood:",object@logLike,"\n")
    cat("BIC:",object@BIC,"\n")
    cat("ICL:",object@ICL,"\n")
    
})
