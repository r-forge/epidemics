\encoding{UTF-8}
\name{epidemics}
\alias{epidemics}
\title{Genetic simulation of epidemics}
\description{
  These functions are under development. Please email the author before
  using them for published results.
}
\usage{
epidemics(n.sample, duration, t.sample = NULL, seq.length = 10000, 
    mut.rate = 1e-05, n.pop = 1, connectivity = NULL, p.disp = 0.1, 
    pop.size = 1e+05, beta, n.ini.inf = 10, t.infectious = 1, 
    t.recover = 2, plot = TRUE, items = c("nsus", "ninf", "nrec"), 
    col = c("blue", "red", grey(0.3)), lty = c(2, 1, 3), pch = c(20, 
        15, 1), file.sizes = "out-popsize.txt", file.sample = "out-sample.txt")
}
\arguments{

}
\value{

}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## SIMPLE EPIDEMIC ##
metapop <- setMetaPop(1, 1e6)
x <- epidemics(n.sample=30, beta=1.1, duration=30, meta=metapop)


## SEVERAL POPULATIONS - 1 DIM STEPPING STONE ##
metapop <- setMetaPop(n.pop = 20, metapop.size = 20e5,
                args.spatial = list(setting="oneDimSS",p.disp = 1e-4))

## show metatpop
plot(metapop, max.cir=0.1, arr=FALSE, annot=FALSE, max.lwd=1)
sampTime <- rep(c(0,20,50,100,150), each=10)

## run simulations
x <- epidemics(n.sample=50, beta=1.1, duration=200, meta=metapop, t.sample=sampTime)

## examine and plot result
x
plot(x)
plot(x,  type="fan", y=fastme.ols)



## MORE COMPLEX SETTINGS ##
## define metapop
metapop <- setMetaPop(n.pop = 100, metapop.size = 100e4,
                args.spatial = list(setting="Delaun",p.disp = 5e-4,connect="rgam",shape=3,rate=1))

plot(metapop, max.cir=0.15, max.lwd=2, arr=FALSE, annot=FALSE)

## run simulations
x <- epidemics(n.sample=50, beta=2, duration=50,meta=metapop,t.sample=0, mut=3e-4)

## plot results
plot(x)
plot(x,  type="fan")

}
