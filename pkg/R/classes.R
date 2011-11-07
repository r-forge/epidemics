########################
## print for 'isolates'
########################
print.isolates <- function(x, ...){
    n <- length(x$gen)
    cat(paste("\n-- ", n, "isolates --"))
    cat("\n- genotypes ($gen) -\n")
    print(x$gen)
    cat("\n- populations ($pop) -\n")
    print(x$pop)
    return(invisible())
}





################
## distIsolates
################
distIsolates <- function(x, add.root = TRUE, res.type = c("dist", "matrix")){
    if (!inherits(x, "isolates"))
        stop("x is not a isolates object")
    res.type <- match.arg(res.type)
    x <- x$gen
    if (add.root) {
        x <- c(list(character(0)), x)
    }
    n <- length(x)
    f1 <- function(a, b){
        return(sum(!union(unlist(a), unlist(b)) %in% intersect(unlist(a),unlist(b))))
    }
    res <- matrix(0, ncol = n, nrow = n)
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            res[i, j] <- f1(x[[i]], x[[j]])
        }
    }
    res <- res + t(res)
    if (res.type == "dist") {
        res <- as.dist(res)
    }
    return(res)
} # end distIsolates






#################
## plot.isolates
#################
plot.isolates <- function(x, y=NULL, ..., plot=TRUE, show.pop=TRUE, col.pal=rainbow, add.root=TRUE){
    ## CHECKS ##
    if(!inherits(x, "isolates")) stop("x is not an 'isolates' object")
    if(!require(ape)) stop("the ape package is needed for phylogenies")
    if(is.null(y)) y <- fastme.ols
    if(!is.function(y)) stop("y argument is not a function")


    ## GET PHYLOGENY ##
    res <- y(distIsolates(x, add.root=add.root, res.type=c("dist")))
    if(!inherits(res, "phylo")) warning("the produced phylogeny is not a 'phylo' object")
    res <- root(res,1)

    ## PLOT RESULTS ##
    if(plot){
        if(show.pop){
            npop <- length(levels(x$pop))
            myCol <- col.pal(npop)[as.integer(x$pop)]
            plot(res, tip.col=myCol, ...)
        } else {
            plot(res, ...)
        }
    }


    ## RETURN TREE ##
    return(invisible(res))
}
