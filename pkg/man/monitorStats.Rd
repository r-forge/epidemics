\encoding{UTF-8}
\name{monitor.epidemics}
\alias{monitor.epidemics}
\title{Monitoring genetic summary statistics during disease epidemics}
\description{
  These functions are under development. Please email the author before
  using them for published results.
}
\usage{
monitor.epidemics(n.sample, duration, seq.length=1e4, mut.rate=1e-5, n.pop=1, connectivity=NULL, p.disp=0.1,
                  pop.size=1e5,  beta, n.ini.inf=10, t.infectious=1, t.recover=2, min.samp.size=100, plot=TRUE
                  items=c("nbSnps","Hs","meanNbSnps","varNbSnps","meanPairwiseDist","varPairwiseDist","meanPairwiseDistStd","varPairwiseDistStd","Fst"),
                  file.sizes="out-popsize.txt", file.sumstat="out-sumstat.txt")
}
\arguments{

}
\value{

}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## SIMPLE EXAMPLE ##
metapop <- setMetaPop(1, 1e6)
monitor.epidemics(300, beta=1.1, duration=100, metaPop=metapop)


## MULTIPLE POPULATIONS - 2D LATTICE ##
metapop <- setMetaPop(n.pop = 20, metapop.size = 20e5,
                args.spatial = list(p.disp = 1e-4))
plot(metapop)

## run simulation
res <- monitor.epidemics(300, beta=1.1, duration=100, min.samp.size=1e3, metaPop=metapop)

## again - keep only certain plots
res <- monitor.epidemics(300, beta=1.1, duration=100, min.samp.size=1e3,
                   metaPop=metapop,item=c("meanPairwiseDistStd", "varPairwiseDistStd", "Fst"))

}
