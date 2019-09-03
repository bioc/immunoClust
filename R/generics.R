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

## label
label <- function(object,...) UseMethod("label")



##
## subset: already S3 in base
setGeneric("subset")

transformParams <- function(object, ... ) UseMethod("transformParams")


## generics only for immunoMeta

## events
events <- function(object,...) UseMethod("events")

prop <- function(object, ...) UseMethod("prop")
"prop<-" <- function(object, ..., value) UseMethod("prop<-")
desc <- function(object, ...) UseMethod("desc")
descFull <- function(object, ...) UseMethod("descFull")
"desc<-" <- function(object, ..., value) UseMethod("desc<-")

level <- function(object, ...) UseMethod("level")
findLevel <- function(object, ...) UseMethod("findLevel")
"level<-" <- function(object, ..., value) UseMethod("level<-")
clusters <- function(object, ...) UseMethod("clusters")
classified <- function(object, ...) UseMethod("classified")
unclassified <- function(object, ...) UseMethod("unclassified")


"addLevel<-" <- function(object, ..., value) UseMethod("addLevel<-")
"move<-" <- function(object, ..., value) UseMethod("move<-")
"parent<-" <- function(object, ..., value) UseMethod("parent<-")
"remove<-" <- function(object, ..., value) UseMethod("remove<-")

"transfer<-" <- function(x, value) UseMethod("transfer<-")
finalize <- function(x,remove.empty=FALSE, depth=-1) UseMethod("finalize")
#makeModel <- function(x,remove.empty=TRUE, depth=-1) UseMethod("makeModel")

## TODO transformParams => transform??
# transform: already S4-method in flowCore
#setGeneric("transform")




