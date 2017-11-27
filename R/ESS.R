#' @name ESS
#' @title ESS for two-player games with a maximum of three strategies
#' @description Computes Evolutionary Stable Strategies of a game with 
#' two players and a maximum of three strategies.
#' @aliases ESS
#' @export ESS
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param A Numeric matrix of size 2x2 or 3x3 representing the number of
#'  strategies of a symmetric matrix game.
#' @param strategies String vector of length n that names all strategies
#'  whereas n represents the number of strategies.
#' @param floats Logical value that handles number representation. If set to
#'  \code{TRUE}, floating-point arithmetic will be used, otherwise fractions.
#'  Default is \code{TRUE}.
#' @return Numeric matrix. Each row represents an ESS.
#' @references Smith, J. M. and Price, G. R. (1973)
#' "The logic of animal conflict", 
#' Nature 246, pp. 15--18.
#' @examples
#' ESS(matrix(c(-1, 4, 0, 2), 2, byrow=TRUE), c("Hawk", "Dove"), FALSE)
#' ESS(matrix(c(1, 2, 0, 0, 1, 2, 2, 0, 1), 3, byrow=TRUE))

ESS <- function(A, strategies = c(), floats = TRUE) {
  if(!is.matrix(A) || !is.numeric(A)) {
    stop("A must be a numeric matrix.")
  }
  else if(nrow(A) != ncol(A)) {
    stop("A must be of size 2x2 or 3x3.")
  }
  else if(nrow(A) != 2 && nrow(A) != 3) {
    stop("A must be of size 2x2 or 3x3.")
  }
  else if(!is.null(strategies) && length(strategies) != ncol(A)) {
    stop("A must be of size nxn.")
  }

  ESS <- c()
  A <- round(A, 5)

  # determine ESS for 2x2 case

  if(nrow(A) == 2) {
    a <- A[1, 1]
    b <- A[1, 2]
    c <- A[2, 1]
    d <- A[2, 2]

    # pure strategy ESS
    if(a > c && b > d) {
      ESS <- rbind(ESS, c(1, 0))
    }
    else if(a > c && b < d) {
      ESS <- rbind(ESS, c(1, 0))
      ESS <- rbind(ESS, c(0, 1))
    }
    else if(a < c && b < d) {
      ESS <- rbind(ESS, c(0, 1))
    }
    # mixed strategy ESS
    else if(a < c && b > d) {
      ESS <- rbind(ESS, c((b - d), (c - a)) / (b + c - a - d))
    }
  }

  # determine ESS for 3x3 case

  # transform matrix in order to have zeros on the diagonal
  if(nrow(A) == 3) {
    for(i in 1:3) {
      if(diag(A)[i] != 0) {
        A[,i] <- A[,i] - diag(A)[i]
      }
    }

    a <- A[1, 2]
    b <- A[1, 3]
    c <- A[2, 1]
    d <- A[2, 3]
    e <- A[3, 1]
    f <- A[3, 2]

    # pure strategy ESS
    if(c < 0 && e < 0) {
      ESS <- rbind(ESS, c(1, 0, 0))
    }
    if(a < 0 && f < 0) {
      ESS <- rbind(ESS, c(0, 1, 0))
    }
    if(b < 0 && d < 0) {
      ESS <- rbind(ESS, c(0, 0, 1))
    }

    # two strategies mixed ESS
    alpha <- a * d + b * f - d * f
    beta <- b * c + d * e - b * e
    gamma <- a * e + c * f - a * c

    if(a > 0 && c > 0 && gamma < 0) {
      ESS <- rbind(ESS, c(a, c, 0) / (a + c))
    }

    if(b > 0 && e > 0 && beta < 0) {
      ESS <- rbind(ESS, c(b, 0, e) / (b + e))
    }

    if(d > 0 && f > 0 && alpha < 0) {
      ESS <- rbind(ESS, c(0, d, f) / (d + f))
    }

    # internal ESS
    firstCond <- (a + c) > 0 && (b + e) > 0 && (d + f) > 0

    if(firstCond) {
      secondCond <- abs(sqrt(b + e) - sqrt(d + f)) < sqrt(a + c) &&
        sqrt(a + c) < (sqrt(b + e) + sqrt(d + f))

      if(secondCond){
        ESS <- rbind(ESS, c(alpha, beta, gamma) / (alpha + beta + gamma))
      }
    }
  }

  # change number representation
  if(!is.null(ESS)) {
    if(!floats) {
      ESS <- MASS::fractions(ESS)
    }

    if(!is.null(strategies)) {
      colnames(ESS) <- strategies
    }
  }

  return(ESS)
}
