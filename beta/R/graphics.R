
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
plot.metaPopInfo <- function(x, y=NULL, ..., col="blue", max.lwd=10, max.cir=0.5, arr=TRUE, annot=TRUE,
                             network.front=TRUE, no.margin=FALSE){
    ## CHECK OBJECT ##
    .check.metaPopInfo(x, stopOnError=FALSE)
    if(x$n.pop==1) return(invisible())

    xy <- x$xy
    cn <- x$cn

    col <- rep(col, length=x$n.pop)

    ## PLOT STUFF ##
    ## empty plot
    rx <- abs(diff(range(xy[,1])))
    ry <- abs(diff(range(xy[,2])))
    xlim <- range(xy[,1]) + c(-0.1, 0.1)*rx
    ylim <- range(xy[,2]) + c(-0.1, 0.1)*ry
    if(no.margin){
        omar <- par("mar")
        on.exit(par(omar))
        par(mar=rep(.1,4))
        plot(xy, , xlab="",ylab="", type="n", xlim=xlim, ylim=ylim, xaxt="n", yaxt="n")

    }else{
        plot(xy, , xlab="x",ylab="y", type="n", xlim=xlim, ylim=ylim)
    }

    ## add pop sizes
    if(network.front)
    symbols(xy, circ=sqrt(x$pop.sizes), xlab="x",ylab="y", inche=max.cir, bg=col, add=TRUE,...)

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

    if(!network.front)
    symbols(xy, circ=sqrt(x$pop.sizes), xlab="x",ylab="y", inche=max.cir, bg=col, add=TRUE,...)

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
## epiMap
##############
epiMap <- function(dat, metapop, groups=c("nsus","ninf","nrec"), sumstat=NULL, max.lwd=3, max.cir=0.2, arr=FALSE, annot=FALSE,
                      out.dir=NULL, ask=TRUE, plot=TRUE, interval=0.05, title="Epidemics simulation",...){
    ## CHECK/PROCESS ARGUMENTS ##
    ## METAPOP PARAMETERS
    .check.metaPopInfo(metapop)

    if(!is.null(out.dir) && !require(animation)) stop("animation package is required - install it or set 'animate' to FALSE.")
    addStat <- !is.null(sumstat)
    if(!all(groups %in% colnames(dat))) stop("some of the requested groups are not in dat")
    if(addStat && !sumstat %in% colnames(dat)) stop("requested summary statistics is not in dat")

    ## GET PARAMETERS TO PLOT ##
    ## xy coords and popsizes ##
    xy <- metapop$xy
    ##temp <- dat[,c("nsus","ninf","nrec")]
    temp <- dat[dat$patch>0,c("nsus","nrec")]
    temp <- prop.table(as.matrix(temp),1)
    inf <- dat$ninf[dat$patch>0]/max(dat$ninf[dat$patch>0])
    if(addStat){
        sumstat.name <- sumstat
        sumstat <- dat[dat$patch>0,sumstat]
        x <- cbind.data.frame(step=dat$step[dat$patch>0], patch=dat$patch[dat$patch>0], inf=inf, temp, sumstat)
        colnames(x)[ncol(x)] <- sumstat.name
    } else {
        x <- cbind.data.frame(step=dat$step[dat$patch>0], patch=dat$patch[dat$patch>0], inf=inf, temp)
        sumstat.name <- NULL
    }
    x.bystep <- split(x, x$step)

    ## define color palette ##
    myPal <- colorRampPalette(c("white","red"))(100)

    ## general dynamics ##
    globalDat <- dat[dat$patch==0,c("step",groups,sumstat.name)]


    ## MAKE SERIES OF PLOTS ##
    par(ask=ask)
    layout(matrix(c(1,2),ncol=1), heights=c(.65,.35))
    if(!is.null(out.dir)) dir.create(out.dir)
    myCol <- c("blue","red","green")
    names(myCol) <- c("nsus","ninf","nrec")
    myCol <- myCol[groups]

    ## OUTPUT TO SCREEN ##
    if(plot){
        for(i in 1:length(x.bystep)){
            ## SYMBOLS ##
            myCirc <- rgb(rep(0,nrow(x.bystep[[i]])), x.bystep[[i]][,"nrec"], x.bystep[[i]][,"nsus"])
            ##myCol <- rgb(1,0,0,x.bystep[[i]][,"inf"])
            temp <- floor(100*x.bystep[[i]][,"inf"])
            temp[temp<1] <- 1
            cirCol <- myPal[temp]
            plot(metapop, col=cirCol, max.lwd=max.lwd, max.cir=max.cir, arr=arr, annot=annot, fg=myCirc, lwd=2,
                 network.front=FALSE, no.margin=TRUE)

            ## GLOBAL DYNAMICS ##
            par(mar=c(4,3,.1,3))
            matplot(globalDat[,"step"], globalDat[,groups], type="n",ylab="Number of individuals", xlab="Time")
            matplot(globalDat[1:i,"step"], globalDat[1:i,groups], col=myCol,lty=1, type="l",lwd=2, add=TRUE)

            ## SUMMARY STATISTICS ##
            if(addStat){
                par(new=TRUE)
                plot(globalDat[,"step"], globalDat[,sumstat.name], xaxt="n",yaxt="n",xlab="",ylab="",type="n")
                points(globalDat[1:i,"step"], globalDat[1:i,sumstat.name], lty=1, type="b",lwd=2)
                axis(side=4)
            }

            ## WAIT IF NEEDED ##
            if(!ask) Sys.sleep(interval)
        }
    }


    ## OUTPUT TO FILE ##
    if(!is.null(out.dir)){

        ## SET ANIMATION PARAMETERS ##
        dir.create(out.dir)
        oani <- ani.options()
        on.exit(ani.options(oani))
        ani.options(nmax=length(x.bystep),outdir=out.dir,interval=0.1,title=title)
        ani.start()

        max.cir <- max.cir/2.5

        for(i in 1:length(x.bystep)){

            layout(matrix(c(1,2),ncol=1), heights=c(.65,.35))

            ## SYMBOLS ##
            myCirc <- rgb(rep(0,nrow(x.bystep[[i]])), x.bystep[[i]][,"nrec"], x.bystep[[i]][,"nsus"])
            ##myCol <- rgb(1,0,0,x.bystep[[i]][,"inf"])
            temp <- floor(100*x.bystep[[i]][,"inf"])
            temp[temp<1] <- 1
            cirCol <- myPal[temp]
            plot(metapop, col=cirCol, max.lwd=max.lwd, max.cir=max.cir, arr=arr, annot=annot, fg=myCirc, lwd=2,
                 network.front=FALSE, no.margin=TRUE)

            ## GLOBAL DYNAMICS ##
            par(mar=c(4,3,.1,3))
            matplot(globalDat[,"step"], globalDat[,groups], type="n",ylab="Number of individuals", xlab="Time")
            matplot(globalDat[1:i,"step"], globalDat[1:i,groups], col=myCol,lty=1, type="l",lwd=2, add=TRUE)

            ## SUMMARY STATISTICS ##
            if(addStat){
                par(new=TRUE)
                plot(globalDat[,"step"], globalDat[,sumstat.name], xaxt="n",yaxt="n",xlab="",ylab="",type="n")
                points(globalDat[1:i,"step"], globalDat[1:i,sumstat.name], lty=1, type="b",lwd=2)
                axis(side=4)
            }


        }

        ani.stop()

    }

}
