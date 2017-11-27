#' @name BNN
#' @title Brown-von Neumann-Nash dynamic
#' @description Brown-von Neumann-Nash replicator dynamic as a type of
#'  evolutionary dynamics.
#' @aliases BNN
#' @export BNN
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param time Regular sequence that represents the time sequence under which
#'  simulation takes place.
#' @param state Numeric vector that represents the initial state.
#' @param parameters Numeric vector that represents parameters needed by the
#'  dynamic.
#' @return Numeric list. Each component represents the rate of change depending on
#'  the dynamic.
#' @references Brown, G. W. and von Neumann, J. (1950) 
#' "Solutions of games by differential equations", In:
#' Kuhn, Harold William and Tucker, Albert William (Eds.) 
#' "Contributions to the Theory of Games I", 
#' Princeton University Press, pp. 73--79.
#' @examples
#' dynamic <- BNN
#' A <- matrix(c(0, -2, 1, 1, 0, -2, -2, 1, 0), 3, byrow=TRUE)
#' state <- matrix(c(0.4, 0.3, 0.3), 1, 3, byrow=TRUE)
#' phaseDiagram3S(A, dynamic, NULL, state, FALSE, FALSE)

BNN <- function(time, state, parameters) {
  a <- parameters
  states <- sqrt(length(a))
  A <- matrix(a, states, byrow = TRUE)
  A <- t(A)

  dX <- c()

  for(i in 1:states) {
    dX[i] <- sum(state * A[i, ])
  }

  avgFitness <- sum(dX * state)

  maxVals <- c()

  for(i in 1:states) {
    maxVals[i] <- max(0, dX[i] - avgFitness)
  }

  for(i in 1:states) {
    dX[i] <- maxVals[i] - state[i] * (sum(maxVals))
  }

  return(list(dX))
}
