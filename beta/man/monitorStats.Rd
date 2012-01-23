\encoding{UTF-8}
\name{monitor.epidemics}
\alias{monitor.epidemics}
\title{Monitoring genetic summary statistics during disease epidemics}
\description{
  These functions are under development. Please email the author before
  using them for published results.
}
\usage{
monitor.epidemics(n.sample, duration, beta, metaPopInfo, seq.length = 10000, 
    mut.rate = 1e-05, n.ini.inf = 10, t.infectious = 1, t.recover = 2, 
    min.samp.size = 100, plot = TRUE, items = c("nbSnps", "Hs", 
        "meanNbSnps", "varNbSnps", "meanPairwiseDist", "varPairwiseDist", 
        "meanPairwiseDistStd", "varPairwiseDistStd", "Fst"), 
    file.sizes = "out-popsize.txt", file.sumstat = "out-sumstat.txt")
}
\arguments{
  \item{n.sample}{the number of samples required.}
  \item{duration}{the duration of the simulation.}
  \item{beta}{the transmission rate.}
  \item{metaPopInfo}{a \code{metaPopInfo} object containing information
    on the infected metapopulation; see \code{\link{setMetaPop}}.}
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
  \item{min.samp.size}{the minimum number of pathogens in the
    metapopulations required to compute summary statistics at a given
    time step.}
  \item{plot}{a logical indicating whether plots should be created at
    the end of the simulation.}
  \item{items}{a vector of character strings indicating which data
    should be plotted.}
  \item{file.sizes}{a character string indicating the name of the output
    file for population dynamics.}
  \item{file.sumstat}{a character string indicating the name of the output
    file for the summary statistics computed.}
}
\value{
  A list containing two slots:

  - \code{$popdyn}: a \code{data.frame} containing the number of susceptible,
  infected, and recovered hosts over time.

  - \code{$sumstat}: a \code{data.frame} containing computed summary statistics.
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## SIMPLE EXAMPLE ##
metapop <- setMetaPop(1, 1e6)
monitor.epidemics(300, beta=1.1, duration=100, metaPop=metapop, n.ini=10)


## MULTIPLE POPULATIONS - 2D LATTICE ##
metapop <- setMetaPop(n.pop = 20, metapop.size = 20e5,
                args.spatial = list(p.disp = 1e-4))
par(mfrow=c(1,1))
plot(metapop)

## run simulation
res <- monitor.epidemics(300, beta=2, duration=100,min.samp.size=1e3,metaPop=metapop, n.ini=10)

## check output
names(res)
head(res$popdyn)
head(res$sumstat)

## visualize each plot separately ##
popdyn <- res$popdyn[res$popdyn$patch==0,]
par(mfrow=c(1,1))
myCol <- colorRampPalette(c("blue","red"))(ncol(res$sumstat)-1)
for(i in 2:ncol(res$sumstat)){
plot(popdyn$step, popdyn$ninf, type="l",lwd=5, col=grey(.7),yaxt="n",ylab="", xlab="")
par(new=TRUE)
plot(res$sumstat$step,res$sumstat[,i],col=myCol[i-1],type="b",lwd=2, pch=20, main=names(res$sumstat)[i])
}

}
