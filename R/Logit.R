#' @name Logit
#' @title Logit dynamic
#' @description Logit dynamic as a type of evolutionary dynamics.
#' @aliases Logit
#' @export Logit
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param time Regular sequence that represents the time sequence under which
#'  simulation takes place.
#' @param state Numeric vector that represents the initial state.
#' @param parameters Numeric vector that represents parameters needed by the
#'  dynamic.
#' @return Numeric list. Each component represents the rate of change depending on
#'  the dynamic.
#' @references Fudenberg, D. and Levine, D. K. (1998)
#' "The Theory of Learning in Games", MIT Press.
#' @examples
#' dynamic <- Logit
#' A <- matrix(c(0, -2, 1, 1, 0, -2, -2, 1, 0), 3, byrow=TRUE)
#' state <- matrix(c(0.4, 0.3, 0.3), 1, 3, byrow=TRUE)
#' eta <- 0.1
#' phaseDiagram3S(A, dynamic, eta, state, FALSE, FALSE)

Logit <- function(time, state, parameters) {
  eta <- parameters[length(parameters)]
  
  a <- parameters[-length(parameters)]
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
    stop("Due to internal restrictions of R, please choose a greater value of 
         eta.")
  }
  
  for(i in 1:states) {
    dX[i] <- (exp(eta^(-1) * dX[i])) / etaVals - state[i]
  }
  
  return(list(dX))
}
