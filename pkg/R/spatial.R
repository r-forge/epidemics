#####################
## Function chooseCN
#####################
chooseCn <- function(xy,ask=TRUE, type=NULL, result.type="nb", d1=NULL, d2=NULL, k=NULL,
                     a=NULL, dmin=NULL, plot.nb=TRUE, edit.nb=FALSE){

  if(is.data.frame(xy)) xy <- as.matrix(xy)
  if(ncol(xy) != 2) stop("xy does not have two columns.")
  if(any(is.na(xy))) stop("NA entries in xy.")
  result.type <- tolower(result.type)
   if(is.null(type) & !ask) stop("Non-interactive mode but no graph chosen; please provide a value for 'type' argument.")

  if(!require(spdep, quiet=TRUE)) stop("spdep library is required.")

  res <- list()

  if(!is.null(d2)){
      if(d2=="dmin"){
          tempmat <- as.matrix(dist(xy))
          d2min <- max(apply(tempmat, 1, function(r) min(r[r>1e-12])))
          d2min <- d2min * 1.0001 # to avoid exact number problem
          d2 <- d2min
      } else if(d2=="dmax"){
          d2max <- max(dist(xy))
          d2max <- d2max * 1.0001 # to avoid exact number problem
          d2 <- d2max
      }
  } # end handle d2

  d1.first <- d1
  d2.first <- d2
  k.first <- k

  ## handle type argument
  if(!is.null(type)){
      type <- as.integer(type)
      if(type < 1 |type > 7) stop("type must be between 1 and 7")
      ask <- FALSE
  }

  ## check for uniqueness of coordinates
  if(any(xyTable(xy)$number>1)){ # if duplicate coords
      DUPLICATE.XY <- TRUE
  } else {
      DUPLICATE.XY <- FALSE
  }


  ## if(is.null(type) & !ask) { type <- 1 }

  ### begin large while ###
  chooseAgain <- TRUE
  while(chooseAgain){
    # re-initialisation of some variables
    d1 <- d1.first
    d2 <- d2.first
    k <- k.first

  ## read type from console
    if(ask){
      temp <- TRUE
      while(temp){
        cat("\nChoose a connection network:\n")
        cat("\t Delaunay triangulation (type 1)\n")
        cat("\t Gabriel graph (type 2)\n")
        cat("\t Relative neighbours (type 3)\n")
        cat("\t Minimum spanning tree (type 4)\n")
        cat("\t Neighbourhood by distance (type 5)\n")
        cat("\t K nearest neighbours (type 6)\n")
        cat("\t Inverse distances (type 7)\n")
        cat("Answer: ")

        type <- as.integer(readLines(n = 1))
        temp <- type < 1 |type > 7
        if(temp) cat("\nWrong answer\n")

        if(type %in% 1:4 & DUPLICATE.XY){
            cat("\n\n== PROBLEM DETECTED ==")
            cat("\nDuplicate locations detected\nPlease choose another graph (5-7) or add random noise to locations (see ?jitter).\n")
            temp <- TRUE
        }

      } # end while
    }
    ##

    ## warning about duplicate xy coords
    if(type %in% 1:4 & DUPLICATE.XY){
        stop("Duplicate locations detected and incompatible with graph type 1-4.\nPlease choose another graph (5-7) or add random noise to locations (see ?jitter).")
    }

    ## graph types
    ## type 1: Delaunay
    if(type==1){
      if(!require(tripack, quiet=TRUE)) stop("tripack library is required.")
      cn <- tri2nb(xy)
    }

    # type 2: Gabriel
    if(type==2){
      cn <- gabrielneigh(xy)
      cn <- graph2nb(cn, sym=TRUE)
    }

    ## type 3: Relative neighbours
    if(type==3){
      cn <- relativeneigh(xy)
      cn <- graph2nb(cn, sym=TRUE)
    }

    ## type 4: Minimum spanning tree
    if(type==4){
      if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
      cn <- ade4::mstree(dist(xy)) # there is also a spdep::mstree
      cn <- neig2nb(cn)
    }

    ## type 5: Neighbourhood by distance
    if(type==5){
      if(is.null(d1) |is.null(d2)){
        tempmat <- as.matrix(dist(xy))
        d2min <- max(apply(tempmat, 1, function(r) min(r[r>1e-12])))
        d2min <- d2min * 1.0001 # to avoid exact number problem
        d2max <- max(dist(xy))
        d2max <- d2max * 1.0001 # to avoid exact number problem
        dig <- options("digits")
        options("digits=5")
        cat("\n Enter minimum distance: ")
        d1 <- as.numeric(readLines(n = 1))
        cat("\n Enter maximum distance \n(dmin=", d2min, ", dmax=", d2max, "): ")
        d2 <- readLines(n = 1)
        ## handle character
        if(d2=="dmin") {
            d2 <- d2min
        } else if(d2=="dmax") {
            d2 <- d2max
        } else {
            d2 <- as.numeric(d2)
        }
        ## restore initial digit option
        options(dig)
      }
    # avoid that a point is its neighbour
      dmin <- mean(dist(xy))/100000
      if(d1<dmin) d1 <- dmin
      if(d2<d1) stop("d2 < d1")
      cn <- dnearneigh(x=xy, d1=d1, d2=d2)
    }

    ## type 6: K nearests
    if(type==6){
      if(is.null(k)) {
        cat("\n Enter the number of neighbours: ")
        k <- as.numeric(readLines(n = 1))
      }
      cn <- knearneigh(x=xy, k=k)
      cn <- knn2nb(cn, sym=TRUE)
    }

    ## type 7: inverse distances
    if(type==7){
        if(is.null(a)) {
            cat("\n Enter the exponent: ")
            a <- as.numeric(readLines(n = 1))
        }
        cn <- as.matrix(dist(xy))
        if(is.null(dmin)) {
            cat("\n Enter the minimum distance \n(range = 0 -", max(cn),"): ")
            dmin <- as.numeric(readLines(n = 1))
        }
        if(a<1) { a <- 1 }
        thres <- mean(cn)/1e8
        if(dmin > thres) dmin <- thres
        cn[cn < dmin] <- dmin
        cn <- 1/(cn^a)
        diag(cn) <- 0
        cn <- prop.table(cn,1)
        plot.nb <- FALSE
        edit.nb <- FALSE
        result.type <- "listw"
    } # end type 7

## end graph types

    if(ask & plot.nb) {
      plot(cn,xy)
      cat("\nKeep this graph (y/n)? ")
    ans <- tolower(readLines(n=1))
      if(ans=="n") {chooseAgain <- TRUE} else {chooseAgain <- FALSE}
    }
    else if(plot.nb){
      plot(cn,xy)
      chooseAgain <- FALSE
    }
  else {chooseAgain <- FALSE}

  }
### end large while

  if(edit.nb) {cn <- edit(cn,xy)}

  if(result.type == "listw") {
      if(type!=7) {
          cn <- nb2listw(cn, style="W", zero.policy=TRUE)
      } else {
          cn <- mat2listw(cn)
          cn$style <- "W"
      }
  }

  res <- cn

  attr(res,"xy") <- xy

  return(res)

} # end chooseCN






##############
## setMetaPop
##############
setMetaPop <- function(n.pop, metapop.size, args.pop.size=list(), args.spatial=list(), check.fail.error=TRUE){
    ## HANDLE ARGUMENTS ##
    args.pop.size$n.pop <- n.pop
    args.pop.size$metapop.size <- metapop.size
    args.spatial$n.pop <- n.pop


    ## GET RESULT ##
    res <- list(n.pop=n.pop, metapop.size=metapop.size)
    res <- c(res, list(pop.sizes=do.call(setPopSizes, args.pop.size)) , do.call(setSpatialConfig, args.spatial))
    names(res$pop.sizes) <- rownames(res$xy) <- names(res$cn) <- names(res$weights) <- 1:n.pop
    colnames(res$xy) <- c('x','y')
    res$call <- match.call()
    class(res) <- "metaPopInfo"


    ## CHECK AND RETURN RESULT ##
    .check.metaPopInfo(res, stopOnError=check.fail.error)
    return(res)
} # end setMetaPop








###############
## setPopSizes
###############
setPopSizes <- function(n.pop, metapop.size, distrib=c("equal","runif","rpois","rgamma"),
                        lambda=NULL, shape=NULL, rate=NULL){
    ## CHECKS ##
    distrib <- match.arg(distrib)


    ## GET POPULATION SIZES ##
    idx <- sample(1:n.pop, 1) # population used for total size adjustment

    ## EQUAL SIZES
    if(distrib=="equal"){
        res <- rep(floor(metapop.size/n.pop),n.pop)
        res[idx] <- res[idx] + metapop.size-sum(res)
    }

    ## UNIFORMLY DISTRIBUTED
    if(distrib=="runif"){
        res <- runif(n.pop,0,1e5)
        res <- round(metapop.size*res/sum(res))
        res[idx] <- res[idx] + metapop.size-sum(res)
    }

    ## POISSON DISTRIBUTED
    if(distrib=="rpois"){
        if(is.null(lambda)) stop("lambda is needed for Poisson-distributed population sizes")
        res <- rpois(n.pop, lambda=lambda)
        res <- round(metapop.size*res/sum(res))
        res[idx] <- res[idx] + metapop.size-sum(res)
    }

    ## GAMMA DISTRIBUTED
    if(distrib=="rgamma"){
        if(is.null(shape) | is.null(rate)) stop("shape and rate are needed for gamma-distributed population sizes")
        res <- rgamma(n.pop, shape=shape, rate=rate)
        res <- round(metapop.size*res/sum(res))
        res[idx] <- res[idx] + metapop.size-sum(res)
    }


    ## RETURN RESULT ##
    if(sum(res)!=metapop.size) warning("final metapopulation size is not the requested size")
    return(res)

} # end setPopSizes








################
## .findLattice
################
.findLattice <- function(n){
    n.row <- floor(sqrt(n))
    n.col <- ceiling(n/n.row)
    while(!isTRUE(all.equal(n, n.row*n.col)))
    {
        n.row <- n.row-1
        n.col <- ceiling(n/n.row)
    }
    return(c(n.row,n.col))
}






####################
## setSpatialConfig
####################
setSpatialConfig <- function(n.pop, setting=c("lattice","Delaunay","Gabriel","oneDimSS","panmix"),
                             n.row=NULL, link=c("rook", "queen"), connect=c("uniform","rgamma"), p.disp=NULL,
                             shape=NULL, rate=NULL, xy=NULL, cn=NULL){
    ## CHECKS ##
    if(!require(spdep, quiet=TRUE)) stop("spdep library is required.")
    if(!require(tripack, quiet=TRUE)) stop("tripack library is required.")
    setting <- match.arg(setting)
    link <- match.arg(link)
    connect <- match.arg(connect)


    ## NON-CUSTOM SETTINGS ##
    if(is.null(xy)){
        ## LATTICE
        if(setting=="lattice"){
            ## get nb of rows and columns
            if(is.null(n.row)){
                temp <- .findLattice(n.pop)
                n.row <- temp[1]
                n.col <- temp[2]
            } else {
                n.col <- ceiling(n.pop/n.row)
                if(!isTRUE(all.equal(n.pop, n.row*n.col))) stop(paste("lattice of" ,n.pop,"populations not possible with", n.row, "rows."))
            }

            ## get connection network and xy coords
            cn <- cell2nb(n.row, n.col, type=link)
            xy <- matrix(as.integer(unlist(strsplit(attr(cn, "region.id"), ":"))), ncol=2, byrow=TRUE)
        }


        ## DELAUNAY TRIANGULATION
        if(setting=="Delaunay"){
            ## get xy coords
            xy <- matrix(runif(n.pop*2, 0, 10), ncol=2)

            ## get network
            cn <- tri2nb(xy)
        }


        ## GABRIEL GRAPH
        if(setting == "Gabriel"){
            cn <- gabrielneigh(xy)
            cn <- graph2nb(cn, sym=TRUE)
        }


        ## ONE DIMENSIONAL STEPING STONE
        if(setting=="oneDimSS"){
            cn <- cell2nb(1, n.pop, type="rook")
            xy <- matrix(as.integer(unlist(strsplit(attr(cn, "region.id"), ":"))), ncol=2, byrow=TRUE)
            xy <- xy[,2:1]
        }

        ## PANMIXY
        if(setting=="panmix"){
            ## get xy coords
            xy <- matrix(runif(n.pop*2, 0, 10), ncol=2)

            ## get network
            cn <- lapply(1:n.pop, function(i) setdiff(1:n.pop,i))
            names(cn) <- 1:n.pop
            class(cn) <- "nb"
        }

    }


    ## GET NETWORK IF ONLY XY PROVIDED ##
    if(!is.null(xy) & is.null(cn)){
        cn <- chooseCn(xy)
    }


    ## GET WEIGHTS ##
    if(connect=="uniform"){
        if(is.null(p.disp)) stop("p.disp is needed for uniform connectivity")
        weights <- lapply(cn, function(e) rep(p.disp/length(e), length(e)))
    }

    if(connect=="rgamma"){
        if(is.null(shape) | is.null(rate)) stop("shape and rate need to be specified for gamma-distributed connectivity")
        if(is.null(p.disp)) stop("p.disp is needed for gamma-distributed connectivity")
        f1 <- function(n){
            out <- rgamma(n, shape=shape,rate=rate)
            out <- p.disp*out/sum(out)
            return(out)
        }
        weights <- lapply(cn, function(e) f1(length(e)))
    }


    ## RETURN RESULT ##
    res <- list(xy=xy, cn=cn, weights=weights)
    return(res)
} # end setSpatialConfig
