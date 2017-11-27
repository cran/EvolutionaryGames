#' @name Replicator
#' @title Replicator dynamic
#' @description Replicator dynamic as a type of evolutionary dynamics.
#' @aliases Replicator
#' @export Replicator
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param time Regular sequence that represents the time sequence under which
#'  simulation takes place.
#' @param state Numeric vector that represents the initial state.
#' @param parameters Numeric vector that represents parameters needed by the
#'  dynamic.
#' @return Numeric list. Each component represents the rate of change depending on
#'  the dynamic.
#' @references Taylor, P. D. and Jonker, L. B. (1978)
#' "Evolutionary stable strategies and game dynamics",
#' Mathematical Biosciences 40 (1-2), pp. 145--156.
#' @examples
#' dynamic <- Replicator
#' A <- matrix(c(0, -2, 1, 1, 0, -2, -2, 1, 0), 3, byrow=TRUE)
#' state <- matrix(c(0.4, 0.3, 0.3), 1, 3, byrow=TRUE)
#' phaseDiagram3S(A, dynamic, NULL, state, FALSE, FALSE)

Replicator <- function(time, state, parameters) {
  a <- parameters
  states <- sqrt(length(a))
  A <- matrix(a, states, byrow = TRUE)
  A <- t(A)

  dX <- c()

  for(i in 1:states) {
    dX[i] <- sum(state * A[i, ])
  }

  avgFitness <- sum(dX * state)

  for(i in 1:states) {
    dX[i] <- state[i] * (dX[i] - avgFitness)
  }

  return(list(dX))
}
