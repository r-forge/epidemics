
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








####################
## plot.metaPopInfo
####################
plot.metaPopInfo <- function(x, y=NULL, ..., col="blue", max.lwd=10, max.cir=0.5, arr=TRUE, annot=TRUE){
    ## CHECK OBJECT ##
    .check.metaPopInfo(x, stopOnError=FALSE)
    if(x$n.pop==1) return(invisible())

    xy <- x$xy
    cn <- x$cn

    col <- seq(col, length=x$n.pop)

    ## PLOT STUFF ##
    ## empty plot
    rx <- abs(diff(range(xy[,1])))
    ry <- abs(diff(range(xy[,2])))
    xlim <- range(xy[,1]) + c(-0.1, 0.1)*rx
    ylim <- range(xy[,2]) + c(-0.1, 0.1)*ry
    plot(xy, , xlab="x",ylab="y", type="n", xlim=xlim, ylim=ylim)

    ## add pop sizes
    symbols(xy, circ=sqrt(x$pop.sizes), xlab="x",ylab="y", inche=max.cir, bg=col, add=TRUE)

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





##############
## mapPopDyn
##############
mapPopDyn <- function(metapop, popdyn, max.lwd=3, max.cir=0.3, arr=TRUE, annot=FALSE){
    ## CHECK/PROCESS ARGUMENTS ##
    ## METAPOP PARAMETERS
    .check.metaPopInfo(metapop)


    ## GET PARAMETERS TO PLOT ##
    xy <- metapop$xy
    splitStep <- split(popdyn, popdyn$step)

    for(i in 1:length(splitStep)){
        pinf <- 
        myCol <-
        plot(metapop, col=myCol)
    }


}
