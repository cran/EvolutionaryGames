#' @name ILogit
#' @title ILogit dynamic
#' @description Imitative Logit dynamic as a type of evolutionary dynamics.
#' @aliases ILogit
#' @export ILogit
#' @author Jochen Staudacher \email{jochen.staudacher@@hs-kempten.de}
#' @param time Regular sequence that represents the time sequence under which
#'  simulation takes place.
#' @param state Numeric vector that represents the initial state.
#' @param parameters Numeric vector that represents parameters needed by the
#'  dynamic.
#' @return Numeric list. Each component represents the rate of change depending on
#'  the dynamic.
#' @references Weibull, J. W. (1997)
#' "Evolutionary Game Theory", MIT Press.
#' @examples
#' dynamic <- ILogit
#' A <- matrix(c(-1, 0, 0, 0, -1, 0, 0, 0, -1), 3, byrow=TRUE)
#' state <- matrix(c(0.1, 0.2, 0.7, 0.2, 0.7, 0.1, 0.9, 0.05, 0.05), 3, 3, byrow=TRUE)
#' eta <- 0.7
#' phaseDiagram3S(A, dynamic, eta, state, TRUE, FALSE)

ILogit <- function(time, state, parameters) {
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
    etaVals <- sum(etaVals, state[i] * exp(eta^(-1) * dX[i]))
  }
  
  if(is.infinite(etaVals)) {
    stop("Due to internal restrictions of R, please choose a greater value of 
         eta.")
  }
  
  for(i in 1:states) {
    dX[i] <- (state[i] * exp(eta^(-1) * dX[i])) / etaVals - state[i]
  }
  
  return(list(dX))
}
