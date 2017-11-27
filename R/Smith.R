#' @name Smith
#' @title Smith dynamic
#' @description Smith dynamic as a type of evolutionary dynamics.
#' @aliases Smith
#' @export Smith
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param time Regular sequence that represents the time sequence under which
#'  simulation takes place.
#' @param state Numeric vector that represents the initial state.
#' @param parameters Numeric vector that represents parameters needed by the
#'  dynamic.
#' @return Numeric list. Each component represents the rate of change depending on
#'  the dynamic.
#' @references Smith, M. J. (1984)
#' "The Stability of a Dynamic Model of Traffic Assignment -- 
#'  An Application of a Method of Lyapunov",
#'  Transportation Science 18, pp. 245--252.
#' @examples
#' dynamic <- Smith
#' A <- matrix(c(0, -2, 1, 1, 0, -2, -2, 1, 0), 3, byrow=TRUE)
#' state <- matrix(c(0.4, 0.3, 0.3), 1, 3, byrow=TRUE)
#' phaseDiagram3S(A, dynamic, NULL, state, FALSE, FALSE)

Smith <- function(time, state, parameters) {
  a <- parameters
  states <- sqrt(length(a))
  A <- matrix(a, states, byrow = TRUE)
  A <- t(A)

  dX <- c()
  val <- c()

  for(i in 1:states) {
    dX[i] <- sum(state * A[i, ])
  }

  fMax <- 0
  sMax <- 0

  for(i in 1:states) {
    for(j in 1:states) {
      if(i == j) next

      fMax <- fMax + state[j] * max(0, dX[i] - dX[j])
      sMax <- sMax + max(0, dX[j] - dX[i])
    }

    val[i] <- fMax - state[i] * sMax
    fMax <- 0
    sMax <- 0
  }

  return(list(val))
}
