#############
## epidemics
#############
epidemics <- function(n.sample, duration, beta, metaPopInfo, t.sample=NULL,
                      seq.length=1e4, mut.rate=1e-5,
                      n.ini.inf=10, t.infectious=1, t.recover=2,
                      plot=TRUE, items=c("nsus", "ninf", "nrec"),
                      col=c("blue", "red", grey(.3)), lty=c(2,1,3), pch=c(2,20,1),
                      file.sizes="out-popsize.txt", file.sample="out-sample.txt"){

    ## CHECK/PROCESS ARGUMENTS ##
    ## METAPOP PARAMETERS
    .check.metaPopInfo(metaPopInfo)

    ## npop
    n.pop <- as.integer(max(metaPopInfo$n.pop[1],1))

    ## cninfo
    cninfo <- .metaPopInfo2cninfo(metaPopInfo)

    ## pop.size
    pop.size <- as.integer(metaPopInfo$pop.sizes)
    if(any(pop.size<1)) stop("pop.size cannot contain values less than 1")


    ## OTHER PARAMETERS
    ## n.sample
    n.sample <- as.integer(max(n.sample[1],1))

    ## duration
    duration <- as.integer(max(duration[1],1))

    ## t.sample
    if(is.null(t.sample)){
        t.sample <- rep(0L, n.sample) # by default, all sampled at the end
    } else {
        if(any(t.sample<0 | t.sample>duration)) stop("t.sample cannot be negative or exceed duration")
        if(length(t.sample) != n.sample) warning("t.sample will be recycled as its length does not match n.sample")
        t.sample <- as.integer(rep(t.sample, length=n.sample))
    }

    ## seq.length
    seq.length <- as.integer(seq.length[1])

    ## mut.rate
    mut.rate <- as.double(max(mut.rate[1],0))
    if(mut.rate < 1e-14) warning("mutation rate is zero")

    ## beta
    beta <- as.double(beta[1])
    if(beta<0) stop("beta (transmission rate) cannot be less than 0")

    ## n.ini.inf
    n.ini.inf <- as.integer(max(n.ini.inf[1],1))

    ## t.infectious
    t.infectious <- as.integer(max(t.infectious,1))

    ## t.recover
    t.recover <- as.integer(max(t.infectious,t.infectious+1))

    ## call run_epidemics ##
    .C("R_epidemics", seq.length, mut.rate, n.pop, pop.size, beta, n.ini.inf, t.infectious, t.recover, n.sample, t.sample, duration, cninfo$nbnb, cninfo$listnb, cninfo$weights, PACKAGE="epidemics")

    ## PLOT ##
    if(plot){
        dat <- read.table("out-popsize.txt", header=TRUE)
        dat <- cbind.data.frame(dat$patch, dat[,items])
        plot.dat <- dat[dat$patch==0, ] # keep only metapop info
        if(any(apply(plot.dat, 1, function(e) all(e<1)))){
            plot.dat <- plot.dat[1:(min(which(apply(plot.dat, 1, function(e) all(e<1))))-1), ]
        }

        ## call to matplot
        matplot(plot.dat, type="b", xlab="time step", ylab="size (number of individuals)", lty=lty, col=col, pch=pch)
        legend("topright", lty=lty, col=col, pch=pch, legend=items)
    }


    ## SAVE POPULATION DYNAMICS ##
    res <- list(popdyn=dat)


    ## GET SAMPLE ##
    if(file.exists("out-sample.txt")){
        txt <- readLines("out-sample.txt")
        res$sample <- list(gen=NULL, pop=NULL)
        res$sample$gen <- txt[seq(2, by=2, length=length(txt)/2)]
        res$sample$gen <- gsub("[[:blank:]]$", "", res$sample$gen)
        res$sample$gen <- lapply(res$sample$gen, function(e) unlist(strsplit(e, " ")))
        res$sample$pop <- factor(txt[seq(1, by=2, length=length(txt)/2)])
        class(res$sample) <- "isolates"
    } else {
        res$sample <- NULL
    }


    ## RENAME FILES ##
    file.rename("out-popsize.txt", file.sizes)
    if(file.exists("out-sample.txt")) file.rename("out-sample.txt", file.sample)

    ## return result ##
    return(res)
} # end epidemics







#####################
## monitor.epidemics
#####################
monitor.epidemics <- function(n.sample, duration, beta, metaPopInfo, seq.length=1e4, mut.rate=1e-5,
                              n.ini.inf=10, t.infectious=1, t.recover=2, min.samp.size=100, plot=TRUE,
                              items=c("nbSnps","Hs","meanNbSnps","varNbSnps","meanPairwiseDist","varPairwiseDist","meanPairwiseDistStd","varPairwiseDistStd","Fst"),
                              file.sizes="out-popsize.txt", file.sumstat="out-sumstat.txt"){

    ## CHECK/PROCESS ARGUMENTS ##
    ## METAPOP PARAMETERS
    .check.metaPopInfo(metaPopInfo)

    ## npop
    n.pop <- as.integer(max(metaPopInfo$n.pop[1],1))

    ## cninfo
    cninfo <- .metaPopInfo2cninfo(metaPopInfo)

    ## pop.size
    pop.size <- as.integer(metaPopInfo$pop.sizes)
    if(any(pop.size<1)) stop("pop.size cannot contain values less than 1")


    ## OTHER PARAMETERS

    ## define useless arguments but needed by the C function
    t.sample <- rep(0L, n.sample) # by default, all sampled at the end

    ## useful arguments ##
    ## n.sample
    n.sample = as.integer(max(n.sample[1],1))

    ## duration
    duration = as.integer(max(duration[1],1))

    ## seq.length
    seq.length <- as.integer(seq.length[1])

    ## mut.rate
    mut.rate <- as.double(max(mut.rate[1],0))
    if(mut.rate < 1e-14) warning("mutation rate is zero")

    ## beta
    beta <- as.double(beta[1])
    if(beta<0) stop("beta (transmission rate) cannot be less than 0")

    ## n.ini.inf
    n.ini.inf <- as.integer(max(n.ini.inf[1],1))

    ## t.infectious
    t.infectious <- as.integer(max(t.infectious,1))

    ## t.recover
    t.recover <- as.integer(max(t.infectious,t.infectious+1))

    ## min.samp.size
    min.samp.size <- as.integer(max(min.samp.size,1))[1]


    ## call R_monitor_epidemics ##
    .C("R_monitor_epidemics", seq.length, mut.rate, n.pop, pop.size, beta, n.ini.inf, t.infectious, t.recover, n.sample, t.sample, duration,
       cninfo$nbnb, cninfo$listnb, cninfo$weights, min.samp.size, PACKAGE="epidemics")


    ## RETRIEVE OUTPUT ##
    ## rename files ##
    sumstat <- read.table("out-sumstat.txt", header=TRUE)
    grpsizes <- read.table("out-popsize.txt", header=TRUE)[,1:4]
    grpsizes <- grpsizes[grpsizes$patch==0,]
    file.rename("out-popsize.txt", file.sizes)
    file.rename("out-sumstat.txt", file.sumstat)

    if(any(apply(grpsizes, 1, function(e) all(e<1)))){
        grpsizes <- grpsizes[1:(min(which(apply(grpsizes, 1, function(e) all(e<1))))-1), ]
    }

    grpsizes <- grpsizes[min(sumstat$step):max(sumstat$step),]


    ## PLOT ##
    if(plot){
        if(length(items)>1 & length(items)<5) par(mfrow = c(2,2))
        if(length(items)>5 & length(items)<7) par(mfrow = c(3,2))
        if(length(items)>6) par(mfrow = c(3,3))
        for(item in items){
            par(mar=c(5,4,3,4))
            matplot(grpsizes$step, grpsizes[,-1], type="l", xlab="time step", col=c("blue", "red", grey(.3)), lty=c(2,1,3), pch=c(2,20,1), yaxt="n", ylab="")
            axis(side=4)
            ##mtext(side=4, "size (number of individuals)", line=2)
            par(new=TRUE)
            plot(sumstat$step, sumstat[,item], ylab=item, xlab="", main=item, type="b", lwd=2)
        }
    }

    ## RETURN RESULT ##
    res <- list(popdyn=grpsizes, sumstat=sumstat)
    return(res)
} # end monitor.epidemics




