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
    if (add.root) { # root added at the end
        x <- c(x, list(character(0)))
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
    res$tip.label <- gsub(" ", "", res$tip.label)
    res <- root(res, as.character(length(res$tip.label)))

    ## PLOT RESULTS ##
    if(plot){
        if(show.pop){
            npop <- length(levels(x$pop))
            ## reorder x$pop
            x$pop <- x$pop[as.numeric(res$tip.label)]
            myCol <- col.pal(npop)[as.integer(x$pop)]
            myCol[is.na(myCol)] <- "black"
            plot(res, tip.col=myCol, ...)
        } else {
            plot(res, ...)
        }
    }


    ## RETURN TREE ##
    return(invisible(res))
}







#####################
## print.metaPopInfo
#####################
.check.metaPopInfo <- function(x, stopOnError=TRUE){
    if(stopOnError) {f1 <- stop} else {f1 <- warning}


    ## GENERAL CHECKS ##
    if(!inherits(x, "metaPopInfo")) f1("metaPopInfo object is not of class metaPopInfo")
    if(!is.list(x)) f1("metaPopInfo object is not a list")


    ## CHECK COMPONENT PRESENCE ##
    x.names <- names(x)
    if(!"n.pop" %in% x.names) f1("metaPopInfo object has no 'n.pop' component.")
    if(!"metapop.size" %in% x.names) f1("metaPopInfo object has no 'metapop.size' component.")
    if(!"pop.sizes" %in% x.names) f1("metaPopInfo object has no 'pop.sizes' component.")
    if(!"xy" %in% x.names) f1("metaPopInfo object has no 'xy' component.")
    if(!"cn" %in% x.names) f1("metaPopInfo object has no 'cn' component.")
    if(!"weights" %in% x.names) f1("metaPopInfo object has no 'weights' component.")
    if(!"call" %in% x.names) f1("metaPopInfo object has no 'call' component.")


    ## CHECK LENGTHS ##
    if(x$n.pop > 1){
        temp <- sapply(x, function(e) ifelse(is.matrix(e), nrow(e), length(e)))
        temp <- temp[names(temp) %in% c("pop.sizes","xy","cn","weights")]
        if(!all(temp==x$n.pop)) f1("inconsistent dimensions found in the content of metaPopInfo object.")
    }


    ## CHECK POP SIZES ##
    if(x$metapop.size != sum(x$pop.sizes)) f1("error in metaPopInfo object: \nsum of population sizes differs from total metapopulation size")

    ## RETURN ##
    return(invisible())
}






#####################
## print.metaPopInfo
#####################
print.metaPopInfo <- function(x, ...){
    cat(paste("\n-- metapopulation with", x$n.pop, "populations --"))
    cat("\n- population sizes ($pop.sizes) -\n")
    print(x$pop.sizes)
    cat("\n- xy coordinates ($xy) -\n")
    print(x$xy)
    cat("\n- connection network ($cn) -\n")
    print(x$cn)
    cat("\n- dispersal probabilities ($weights) -\n")
    print(x$weights)
    cat("\n- call ($call) -\n")
    print(x$call)
    return(invisible())
} # end print.metaPopInfo







####################
## plot.metaPopInfo
####################
plot.metaPopInfo <- function(x, y=NULL, ..., max.lwd=10, max.cir=0.5, arr=TRUE, annot=TRUE){
    ## CHECK OBJECT ##
    .check.metaPopInfo(x, stopOnError=FALSE)
    if(x$n.pop==1) return(invisible())

    xy <- x$xy
    cn <- x$cn

    ## PLOT STUFF ##
    ## empty plot
    rx <- abs(diff(range(xy[,1])))
    ry <- abs(diff(range(xy[,2])))
    xlim <- range(xy[,1]) + c(-0.1, 0.1)*rx
    ylim <- range(xy[,2]) + c(-0.1, 0.1)*ry
    plot(xy, , xlab="x",ylab="y", type="n", xlim=xlim, ylim=ylim)

    ## add pop sizes
    symbols(xy, circ=sqrt(x$pop.sizes), xlab="x",ylab="y", inche=max.cir, bg="blue3", add=TRUE)

    ## add arrows
    f1 <- function(i){ # find arrow parameters
        startx <- rep(xy[i,1], length(cn[[i]]))
        starty <- rep(xy[i,2], length(cn[[i]]))
        endx <- xy[cn[[i]], 1]
        endy <- xy[cn[[i]], 2]
        w <- x$weights[[i]]
        res <- data.frame(startx, starty, endx, endy, w)
        return(res)
    }

    temp <- Reduce("rbind.data.frame",lapply(1:length(x$cn), f1))
    arr.w <- max.lwd*temp$w/max(temp$w)
    if(arr) {
        arrows(temp$startx, temp$starty, temp$endx, temp$endy, lwd=arr.w)
    } else {
        segments(temp$startx, temp$starty, temp$endx, temp$endy, lwd=arr.w)
    }

    ## add labels
    if(annot){
        text(xy, lab=1:nrow(xy), col="grey", cex=1.5, font=2)
        text(xy, lab=1:nrow(xy), col="white", cex=1.5)
    }


    ## RETURN RES ##
    res <- list(xy=xy, circ=sqrt(x$pop.sizes), arrows=temp)
    return(invisible(res))
}
