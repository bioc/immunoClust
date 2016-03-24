####
##
## perform asinh transformation on flourescence parameters 
##
####


### trans.ApplyToData
trans.ApplyToData <- function(
x, data, add.param=c(), 
max.decade=attr(x,"trans.decade"), lin.scale=attr(x,"trans.scale") 
) {
    
    ret <- data
    
    if( !is.null(x@label) & length(x@label) < nrow(data) ) {
        ret <- data[1:length(x@label),]
    }
    
    if( is.null(max.decade) ) {
        max.decade = -1
    }
    
    if( is.null(lin.scale) ) {
        lin.scale <- 1.0
    }
    
    param <- attributes(x)$param
    
    if( class(data) == "flowFrame" ) {
        
        mat <- as.matrix(exprs(data))[,c(param, add.param)]
        
        par <- parameters(data)
        range <- range(data)[c(param,add.param)]
        P <- length(param)
        for( i in 1:P ) {
            if( x@trans.a[i] > 0.0 ) {
                a <- x@trans.a[i]
                b <- x@trans.b[i]
                x_max <- max(mat[,i])
                C <- 1.0
                if( max.decade > 0 ) {
                    C <- max.decade / asinh(a * x_max + b)
                }
                mat[,i] <- C * asinh(a * mat[,i] + b )
                j <- match(param[i], colnames(data))
                par$maxRange[[j]] <- C * asinh(a*range[2,i] + b)
                par$minRange[[j]] <- min(mat[,i])
            }
            else {
                mat[,i] <- mat[,i]/lin.scale
                j <- match(param[i], colnames(data))
                par$maxRange[[j]] <- range[2,i]/lin.scale
                par$minRange[[j]] <- range[1,i]/lin.scale
                
            }
            
        }
        
        
        inc <- match(c(param,add.param), par@data[,'name'])
        par@data <- par@data[inc,]
        attr(mat, "ranges") <- NULL
        
        ret <- new("flowFrame", exprs=mat, parameters=par,
                description=description(data))
    }
    else {
        mat <- data[,param]
        P <- length(param)
        for( i in 1:P ) if( x@trans.a[i] > 0.0 ) {
            a <- x@trans.a[i]
            b <- x@trans.b[i]
            x_max <- max(mat[,i])
            C <- 1.0
            if( max.decade > 0 ) {
                C <- max.decade / asinh(a * x_max + b)
            }            
            mat[,i] <- C * asinh(a * mat[,i] + b )
        }
        
        ret[,param] <- mat
    }
    
    ret
}
# trans.ApplyToData

### trans.FitToData
trans.FitToData <- function( 
x, data, B=10, tol=1e-5, certainty=0.3, proc="vsHtransAw"
) {

    inc <- !is.na(x@label)
    y <- .exprs(data, x@parameters)[inc,]
    z <- as.matrix(x@z[inc,])
    
    
    N <- nrow(y)
    P <- ncol(y)
    K <- x@K
    
    obj <- .Call(paste(sep="", "immunoC_", proc), 
                as.integer(N), as.integer(P), as.integer(K),
                as.double(t(y)),as.double(t(z)),
                a=as.double(x@trans.a), b=as.double(x@trans.b),
                as.integer(B), as.double(tol), as.double(certainty))
    
    r <- rbind(obj$a,obj$b)
    colnames(r) <- x@parameters
    
    r
    
}
# trans.FitToData
