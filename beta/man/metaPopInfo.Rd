\encoding{UTF-8}
\name{setMetaPop}
\alias{setMetaPop}
\alias{setPopSizes}
\alias{setSpatialConfig}
\alias{metaPopInfo}
\alias{print.metaPopInfo}
\alias{plot.metaPopInfo}
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

setPopSizes(n.pop, metapop.size, distrib = c("equal", "runif",
    "rpois", "rgamma"), lambda = NULL, shape = NULL, rate = NULL)

setSpatialConfig(n.pop, setting = c("lattice", "Delaunay", "Gabriel", 
    "oneDimSS", "panmix"), n.row = NULL, link = c("rook", "queen"), 
    connect = c("uniform", "rgamma"), p.disp = NULL, shape = NULL, 
    rate = NULL, xy = NULL, cn = NULL)

\method{print}{metaPopInfo}(x, \dots)

\method{plot}{metaPopInfo}(x, y = NULL, \dots, max.lwd = 10, max.cir = 0.5, arr = TRUE, 
    annot = TRUE)

}
\arguments{
  \item{n.pop}{the number of populations in the metapopulation.}
  \item{metapop.size}{the total size of the metapopulation.}
  \item{args.pop.size}{a list of named arguments passed to \code{setPopSizes}.}
  \item{args.spatial}{a list of named arguments passed to \code{setSpatialConfig}.}
  \item{check.fail.error}{a logical indicating if in case of issue
    detected in the output, the function should stop with an error
    message (alternative gives warnings only).}
  \item{distrib}{a character string indicating the distribution of
    population sizes - see details.}
  \item{lambda}{parameter of a Poisson distribution.}
  \item{shape, rate}{parameters of a gamma distribution.}
  \item{setting}{a character string indicating the spatial distribution of
    the populations - see details.}
  \item{n.row}{the number of rows to be used for the lattice.}
  \item{link}{the type of link to be used for the lattice, 'rook' or
    'queen'.}
  \item{connect}{a character string indicating the distribution of the
    dispersal probabilities for each population - see details.}
  \item{p.disp}{the probability of dispersal from any given population.}
  \item{xy}{a matrix with two columns containing spatial (x,y) coordinates.}
  \item{cn}{a connectivity network with class \code{nb} (\code{spdep} package).}
}
\value{
  For \code{setMetaPop}, a \code{metaPopInfo} object.

  For \code{chooseCn}, a connectivity network.

  For \code{plot.metaPopInfo}, a \code{list} containing plot information.

}
\details{
  == DISTRIBUTION OF POPULATION SIZES ==\cr
  \code{distrib} indicates the distribution of population sizes in the
  metapopulation; it can take the values:\cr
  - "equal": all populations have the same size (if possible).\cr
  - "runif": population sizes are drawn from a uniform distribution.\cr
  - "rpois": population sizes are drawn from a Poisson distribution with
  parameter \code{lambda}.\cr
  - "rgamma": population sizes are drawn from a gamma distribution with
  parameters \code{shape} and \code{rate}.\cr
  - - - - - - - - - - - - - - - - - -\cr
  == SPATIAL DISTRIBUTION OF THE POPULATIONS ==\cr
  \code{setting} indicates the spatial distribution of the populations;
  it can take the following values:\cr
  - "lattice": populations are arranged on a lattice, aiming for a
  square lattice by default, unless \code{n.row} is provided.\cr
  - "Delaunay": populations locations are drawn from uniform
  distributions, and links are defined according to Delaunay
  triangulation.\cr
  - "Gabriel": populations locations are drawn from uniform
  distributions, and links are defined according to the graph of
  Gabriel.\cr
  - "oneDimSS": populations are arranged on a one-dimensional stepping stone.\cr
  - "panmix": all populations are connected (i.e., panmixie).\cr
  - - - - - - - - - - - - - - - - - -\cr
  == DISTRIBUTION OF THE DISPERSAL PROBABILITIES ==\cr
  \code{connect} indicates the distribution of the dispersal
  probabilities for any given population. Probabilities are rescaled to
  ensure that the overal dispersal probability equals
  \code{p.disp}.\cr
  Possible values are:\cr
  - "uniform": dispersal probabilities are equal for all connections.\cr
  - "rgamma": dispersal probabilities are drawn from a gamma
  distribution with parameters \code{shape} and \code{rate}.
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
