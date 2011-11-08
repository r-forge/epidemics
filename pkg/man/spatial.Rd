\encoding{UTF-8}
\name{setMetaPop}
\alias{setMetaPop}
\alias{metaPopInfo}
\alias{print.metaPopInfo}
\alias{plot.metaPopInfo}
\alias{chooseCn}
\alias{setPopSizes}
\alias{setSpatialConfig}
% \alias{}
% \alias{}
\title{Generate metapopulation settings for epidemic simulation}
\description{
  These functions are under development. Please email the author before
  using them for published results.
}
\usage{
setMetaPop(n.pop, metapop.size, args.pop.size=list(),
           args.spatial=list(), check.fail.error=TRUE)

setPopSizes(n.pop, metapop.size, distrib=c("equal","runif","rpois","rgamma"),
            lambda=NULL, shape=NULL, rate=NULL)

setSpatialConfig(n.pop, setting=c("lattice","Delaunay","Gabriel","oneDimSS","panmix"),
                 n.row=NULL, link=c("rook", "queen"), connect=c("uniform","rgamma"),
                 p.disp=NULL, shape=NULL, rate=NULL, xy=NULL, cn=NULL)

}
\arguments{

}
\value{

}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## 10 pop in a lattice ##
x <- setMetaPop(n.pop = 10, metapop.size = 1e5, args.spatial = list(p.disp = 0.001))
x
plot(x)


## 36 pop in a square lattice, queen relationship, varied pop sizes and connectivity ##
x <- setMetaPop(n.pop = 36, metapop.size = 1e7,
                args.spatial = list(p.disp = 0.001,connect="rgam",shape=1,rate=2,link="queen"),
                args.pop=list(dist="runif"))
x
plot(x, max.cir=0.3)


## Delaunay ##
x <- setMetaPop(n.pop = 20, metapop.size = 1e5,
                args.spatial = list(setting="Delau",p.disp = 1e-4,connect="rgam",shape=1,rate=1),
                args.pop=list(dist="rpois", lambda=2))
x
plot(x, max.cir=0.4, max.lwd=3)


## one-dim stepping stone, 15 pop ##
x <- setMetaPop(n.pop = 8, metapop.size = 8e5,
                args.spatial = list(setting="oneDimSS",p.disp = 1e-4))
x
plot(x, max.cir=0.15, max.lwd=2)


## big Delaunay setting ##
x <- setMetaPop(n.pop = 300, metapop.size = 300e4,
                args.spatial = list(setting="Delaun",p.disp = 1e-4,connect="rgam",shape=3,rate=1),
                args.pop=list(dist="rpois", lambda=3))
x
plot(x, max.cir=0.15, max.lwd=2, arr=FALSE, annot=FALSE)


}
