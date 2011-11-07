\encoding{UTF-8}
\name{monitor.epidemics}
\alias{monitor.epidemics}
\title{Monitoring genetic summary statistics during disease epidemics}
\description{
  These functions are under development. Please email the author before
  using them for published results.
}
\usage{
monitor.epidemics(n.sample, duration, seq.length = 10000, mut.rate = 1e-05, 
    n.pop = 1, connectivity = NULL, p.disp = 0.1, pop.size = 1e+05, 
    beta, n.ini.inf = 10, t.infectious = 1, t.recover = 2, min.samp.size = 100, 
    plot = TRUE, items = c("nbSnps", "Hs", "meanNbSnps", "varNbSnps", 
        "meanPairwiseDist", "varPairwiseDist", "Fst"), file.sizes = "out-popsize.txt", 
    file.sumstat = "out-sumstat.txt")
}
\arguments{

}
\value{

}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## simple example
monitor.epidemics(300, beta=1.1, duration=100)

## multiple populations - several cycles
res <- monitor.epidemics(300, beta=1.1, n.pop=20, duration=100,min.samp.size=1e4, p.disp=2e-5)

## various plots
res <- monitor.epidemics(300, beta=1.1, n.pop=20, duration=100,min.samp.size=1e4, p.disp=2e-5, item=c("meanPairwiseDist","meanPairwiseDistStd", "varPairwiseDist", "varPairwiseDistStd"))
}
