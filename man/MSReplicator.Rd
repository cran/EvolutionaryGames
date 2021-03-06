% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSReplicator.R
\name{MSReplicator}
\alias{MSReplicator}
\title{Maynard Smith replicator dynamic}
\usage{
MSReplicator(time, state, parameters)
}
\arguments{
\item{time}{Regular sequence that represents the time sequence under which
simulation takes place.}

\item{state}{Numeric vector that represents the initial state.}

\item{parameters}{Numeric vector that represents parameters needed by the
dynamic.}
}
\value{
Numeric list. Each component represents the rate of change depending on
 the dynamic.
}
\description{
Maynard Smith replicator dynamic as a type of evolutionary
 dynamics.
}
\examples{
dynamic <- MSReplicator
A <- matrix(c(0, -2, 1, 1, 0, -2, -2, 1, 0), 3, byrow=TRUE)
state <- matrix(c(0.4, 0.3, 0.3), 1, 3, byrow=TRUE)
phaseDiagram3S(A, dynamic, NULL, state, FALSE, FALSE)
}
\references{
Smith, J. M. (1982) 
"Evolution and the Theory of Games",
Cambridge University Press.
}
\author{
Daniel Gebele \email{dngebele@gmail.com}
}
