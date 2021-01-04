## generic
# nobs already S3 in stats
setGeneric("nobs")

npar <- function(object, ...) UseMethod("npar")
ncls <- function(object, ...) UseMethod("ncls")
## weights already S3 in stats
setGeneric("weights")
## mu
mu <- function(object,...) UseMethod("mu")

## sigma already S3 in stats
setGeneric("sigma")
## parameters already S4 in flowCore
setGeneric("parameters")
## parameters<- already S4 in flowCore
#"parameters<-" <- function(object, ..., value) UseMethod("parameters<-")
#setGeneric("parameters<-")

## label
label <- function(object,...) UseMethod("label")

## posterior
posterior <- function(object, ...) UseMethod("posterior")

##
## subset: already S3 in base
setGeneric("subset")

transformParams <- function(object, ... ) UseMethod("transformParams")


## generics only for immunoMeta S3methods (replacement)

## events
events <- function(object,...) UseMethod("events")

prop <- function(object, ...) UseMethod("prop")
"prop<-" <- function(object, ..., value) UseMethod("prop<-")
setGeneric("prop<-")
desc <- function(object, ...) UseMethod("desc")
descFull <- function(object, ...) UseMethod("descFull")
"desc<-" <- function(object, ..., value) UseMethod("desc<-")
setGeneric("desc<-")

level <- function(object, ...) UseMethod("level")
findLevel <- function(object, ...) UseMethod("findLevel")
"level<-" <- function(object, ..., value) UseMethod("level<-")
setGeneric("level<-")
clusters <- function(object, ...) UseMethod("clusters")
classified <- function(object, ...) UseMethod("classified")
unclassified <- function(object, ...) UseMethod("unclassified")


"addLevel<-" <- function(object, ..., value) UseMethod("addLevel<-")
setGeneric("addLevel<-")
"move<-" <- function(object, ..., value) UseMethod("move<-")
setGeneric("move<-")

"parent<-" <- function(object, ..., value) UseMethod("parent<-")
setGeneric("parent<-")
"remove<-" <- function(object, ..., value) UseMethod("remove<-")
setGeneric("remove<-")

"transfer<-" <- function(object, ..., value) UseMethod("transfer<-")
setGeneric("transfer<-")
finalize <- function(object, ...) UseMethod("finalize")
#makeModel <- function(x,remove.empty=TRUE, depth=-1) UseMethod("makeModel")

nsam <- function(object, ...) UseMethod("nsam")
sam_ncls <- function(object, ...) UseMethod("sam_ncls")
sam_clsWeights <- function(object, ...) UseMethod("sam_clsWeights")
sam_clsEvents <- function(object, ...) UseMethod("sam_clsEvents")
sam_clsMu <- function(object, ...) UseMethod("sam_clsMu")
sam_clsSigma <- function(object, ...) UseMethod("sam_clsSigma")

clusterDist <- function(object, ...) UseMethod("clusterDist")
clusterCoeff <- function(object, ...) UseMethod("clusterCoeff")
clusterProb <- function(object, ...) UseMethod("clusterProb")

## TODO transformParams => transform??
# transform: already S4-method in flowCore
#setGeneric("transform")




