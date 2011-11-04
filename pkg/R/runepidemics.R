#############
## epidemics
#############
epidemics <- function(n.sample, duration, t.sample=NULL,
                      seq.length=1e4, mut.rate=1e-5, n.pop=1, connectivity=NULL, p.disp=0.01,
                      pop.sizes,  beta, n.ini.inf=10, t.infectious=1, t.recover=2){

    ## check/process arguments ##
    ## n.sample
    n.sample = as.integer(max(n.sample[1],1))

    ## duration
    duration = as.intger(max(duration[1],1))

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

    ## npop
    n.pop <- as.integer(max(n.pop[1],1))

    ## connectivity
    if(is.null(connectivity)){
        if(p.disp[1]<0 | p.disp[1]>1) stop("pdisp must be comprised between 0 and 1")
        connectivity <- as.double(diag(p.disp, n.pop))
    } else {
        if(!is.matrix(connectivity)) stop("connectivity must provided as a matrix")
        if(nrow(connectivity) != ncol(connectivity)) stop("connectivity matrix must be square")
        if(any(connectivity)<0) stop("connectivity matrix contains negative values")
        connectivity <- as.double(t(connectivity))
    }

    ## pop.size
    pop.size <- as.integer(rep(pop.size, length=n.pop))
    if(any(pop.size<1)) stop("pop.size cannot contain values less than 1")

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
    .C("epidemics", seq.length, mut.rate, n.pop, pop.sizes, beta, n.ini.inf, t.infectious, t.recover, n.sample, t.sample, duration, connectivity)

    ## return result ##

}
