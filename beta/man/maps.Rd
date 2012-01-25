\encoding{UTF-8}
\name{epiMap}
\alias{epiMap}
\title{Generate maps of population dynamics and genetic diversity}
\description{
  These functions are under development. Please email the author before
  using them for published results.
}
\usage{
epiMap(popdyn, metapop, max.lwd=3, max.cir=0.2, arr=FALSE,
         annot=FALSE, out.dir=NULL, ask=TRUE, plot=TRUE,\dots)
}
\arguments{
 
}
\value{
 

}

\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{

#### ONLY POPULATION DYNAMICS ####
## 15 pop in a one dim stepping stone ##
metapop <- setMetaPop(n.pop = 15, metapop.size = 1e6, args.spatial = list(setting="oneDimSS",p.disp = 0.001))
plot(metapop)

x <- epidemics(duration=50, beta=2, metaPop=metapop,n.samp=200)

epiMap(x$popdyn,metapop,ask=FALSE)


## Delaunay ##
metapop <- setMetaPop(n.pop = 300, metapop.size = 1e6,
                args.spatial = list(setting="Delaunay",p.disp = 1e-2),
                args.pop=list(dist="runif"))

plot(metapop, max.cir=0.1, max.lwd=2, network.front=FALSE,annot=FALSE, arr=FALSE)

x <- epidemics(duration=120, beta=2, metaPop=metapop,n.samp=200)

epiMap(x$popdyn,metapop,plot=FALSE,out.dir="tempfig")



#### WITH SUMMARY STATISTICS ####
metapop <- setMetaPop(n.pop = 100, metapop.size = 1e6,
                     args.spatial = list(setting="Delaunay",p.disp = 1e-3))
par(mfrow=c(1,1))
plot(metapop, max.cir=.2, arr=FALSE, network.front=FALSE, annot=FALSE)

## run simulation
res <- monitor.epidemics(150, beta=2, duration=100,min.samp.size=1e3,metaPop=metapop, n.ini=10)

## map dynamic map - number of SNPs
epiMap(res,metapop, "nbSnps")
epiMap(res,metapop, "Fst")

if(require(animation)){
epiMap(res,metapop, "meanPairwiseDistStd",plot=FALSE, out.dir="anim-temp")
}
}
