\encoding{UTF-8}
\name{epidemics}
\alias{epidemics}
\title{Genetic simulation of epidemics}
\description{
  These functions are under development. Please email the author before
  using them for published results.
}
\usage{
epidemics(n.sample, duration, beta, metaPopInfo, t.sample = NULL, 
    seq.length = 10000, mut.rate = 1e-05, n.ini.inf = 10, t.infectious = 1, 
    t.recover = 2, plot = TRUE, items = c("nsus", "ninf", "nrec"), 
    col = c("blue", "red", grey(0.3)), lty = c(2, 1, 3), pch = c(20, 
        15, 1), file.sizes = "out-popsize.txt", file.sample = "out-sample.txt") 
}
\arguments{
  \item{n.sample}{the number of samples required.}
  \item{duration}{the duration of the simulation.}
  \item{beta}{the transmission rate.}
  \item{metaPopInfo}{a \code{metaPopInfo} object containing information
    on the infected metapopulation; see \code{\link{setMetaPop}}.}
  \item{t.sample}{a vector of sampling dates of length \code{n.sample},
    expressed in number of time steps from the end of the simulations;
    only positive integers are accepted.}
  \item{seq.length}{the length of the pathogenic genome, in number of nucleotides.}
  \item{mut.rate}{the mutation rate of the pathogenic genome, in number
    of mutation per site and per time step.}
  \item{n.ini.inf}{the initial number of (identical) pathogens seeding
    the original infection. Used to avoid stochatistc extinctions during
    the early stages of the epidemic.}
  \item{t.infectious}{the age, in number of time steps, at which
    pathogens start being infectious.}
  \item{t.recover}{the age, in number of time steps, at which
    pathogens stops being infectious, causing the host to move to
    'recovered' state.}
  \item{plot}{a logical indicating whether plots should be created at
    the end of the simulation.}
  \item{items}{a vector of character strings indicating which data
    should be plotted.}
  \item{col,lty,pch}{graphical parameters indicating the color, line
    type, and type of point to be used in the plot.}
  \item{file.sizes}{a character string indicating the name of the output
    file for population dynamics.}
  \item{file.sample}{a character string indicating the name of the output
    file for the sampled isolates.}
}
\value{
  A list containing two slots:

  - \code{$popdyn}: a \code{data.frame} containing the number of susceptible,
  infected, and recovered hosts over time.

  - \code{$sample}: a list of class \code{isolates} containing the sampled isolates.
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\seealso{
  \code{\link{setMetaPop}} to set metapopulation parameters.
}
\examples{
## SIMPLE EPIDEMIC ##
metapop <- setMetaPop(1, 1e6)
x <- epidemics(n.sample=30, beta=1.5, duration=100, meta=metapop, n.ini=10)


## SEVERAL POPULATIONS - 1 DIM STEPPING STONE ##
metapop <- setMetaPop(n.pop = 20, metapop.size = 20e5,
                args.spatial = list(setting="oneDimSS",p.disp = 1e-4))

## show metatpop
plot(metapop, max.cir=0.1, arr=FALSE, annot=FALSE, max.lwd=1)
sampTime <- rep(c(0,20,50,100,150), each=10)

## run simulations
x <- epidemics(n.sample=50, beta=2, duration=200,meta=metapop,t.sample=sampTime, n.ini=10)

## examine and plot result
class(x)
names(x)
lapply(x,head)
plot(x$sample)
plot(x$sample,  type="fan", y=fastme.ols)



## MORE COMPLEX SETTINGS ##
## define metapop
metapop <- setMetaPop(n.pop = 100, metapop.size = 100e4,
                args.spatial = list(setting="Delaun",p.disp = 5e-4,connect="rgam",shape=3,rate=1))

plot(metapop, max.cir=0.15, max.lwd=2, arr=FALSE, annot=FALSE)

## run simulations
x <- epidemics(n.sample=50, beta=2,duration=100,meta=metapop,t.sample=0, mut=3e-4, n.ini=10)

## plot results
plot(x$sample)
plot(x$sample,  type="fan")

}
