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
## simple epidemic
epidemics(n.sample=30, beta=1.1, duration=30, t.sample=15)

## several populations - multiple cycles
epidemics(n.sample=30, beta=1.1, duration=200, n.pop=10,p.disp=2e-5)

}
