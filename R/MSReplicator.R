#' @name MSReplicator
#' @title Maynard Smith replicator dynamic
#' @description Maynard Smith replicator dynamic as a type of evolutionary
#'  dynamics.
#' @aliases MSReplicator
#' @export MSReplicator
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param time Regular sequence that represents the time sequence under which
#'  simulation takes place.
#' @param state Numeric vector that represents the initial state.
#' @param parameters Numeric vector that represents parameters needed by the
#'  dynamic.
#' @return Numeric list. Each component represents the rate of change depending on
#'  the dynamic.
#' @references Smith, J. M. (1982) 
#' "Evolution and the Theory of Games",
#' Cambridge University Press.
#' @examples
#' dynamic <- MSReplicator
#' A <- matrix(c(0, -2, 1, 1, 0, -2, -2, 1, 0), 3, byrow=TRUE)
#' state <- matrix(c(0.4, 0.3, 0.3), 1, 3, byrow=TRUE)
#' phaseDiagram3S(A, dynamic, NULL, state, FALSE, FALSE)

MSReplicator <- function(time, state, parameters) {
  a <- parameters
  states <- sqrt(length(a))
  A <- matrix(a, states, byrow = TRUE)
  A <- t(A)
  minVal <- min(A)

  # rescale payoff matrix
  if(minVal <= 0) {
    A <- A + (1 + abs(minVal))
  }

  dX <- c()

  for(i in 1:states) {
    dX[i] <- sum(state * A[i, ])
  }

  avgFitness <- sum(dX * state)

  for(i in 1:states) {
    dX[i] <- (state[i] * (dX[i] - avgFitness)) / avgFitness
  }

  return(list(dX))
}
