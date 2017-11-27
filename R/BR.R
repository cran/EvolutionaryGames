#' @name BR
#' @title BR dynamic
#' @description Best response dynamic as a type of evolutionary dynamics.
#' @aliases BR
#' @export BR
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param time Regular sequence that represents the time sequence under which
#'  simulation takes place.
#' @param state Numeric vector that represents the initial state.
#' @param parameters Numeric vector that represents parameters needed by the
#'  dynamic.
#' @return Numeric list. Each component represents the rate of change depending on
#'  the dynamic.
#' @references Gilboa, I. and Matsui, A. (1991)
#' "Social Stability and Equilibrium",
#' Econometrica 59, pp. 859--867.
#' @examples
#' dynamic <- BR
#' A <- matrix(c(0, -2, 1, 1, 0, -2, -2, 1, 0), 3, byrow=TRUE)
#' state <- matrix(c(0.4, 0.3, 0.3), 1, 3, byrow=TRUE)
#' phaseDiagram3S(A, dynamic, NULL, state, FALSE, FALSE)

BR <- function(time, state, parameters) {
  eta <- 0.002
  a <- parameters
  states <- sqrt(length(a))
  A <- matrix(a, states, byrow = TRUE)
  A <- t(A)
  
  dX <- c()
  
  for(i in 1:states) {
    dX[i] <- sum(state * A[i, ])
  }
  
  etaVals <- c()
  
  for(i in 1:states) {
    etaVals <- sum(etaVals, exp(eta^(-1) * dX[i]))
  }
  
  if(is.infinite(etaVals)) {
    stop("Due to internal restrictions of R, unable to simulate model.")
  }
  
  for(i in 1:states) {
    dX[i] <- (exp(eta^(-1) * dX[i])) / etaVals - state[i]
  }
  
  return(list(dX))
}
