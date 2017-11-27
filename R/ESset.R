#' @name ESset
#' @title Evolutionarily stable set for two-player games with three strategies
#' @description Computes evolutionarily stable sets of a game with two players 
#' and three strategies.
#' @aliases ESset
#' @export ESset
#' @references Thomas, B. (1985)
#' "On evolutionarily stable sets", 
#' Journal of Mathematical Biology 22, pp. 105--115.
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param A Numeric matrix of size 3x3 representing the number of strategies of
#'  a symmetric matrix game.
#' @param strategies String vector of length 3 that names all strategies.
#' @param floats Logical value that handles number representation. If set to
#'  \code{TRUE}, floating-point arithmetic will be used, otherwise fractions.
#'  Default is \code{TRUE}.
#' @return Numeric matrix. Each row represents the start and end point of a
#'  line (ESset). In addition, a plot of the ESset in the game will be created.
#' @examples
#' # Please note that the computation of evolutionarily stable sets 
#' # is rather time-consuming. 
#' # Depending on your machine you might need to wait more 
#' # than 10 seconds in order to run the following example.
#' \dontrun{
#' A <- matrix(c(-2, 5, 10/9, 0, 5/2, 10/9, -10/9, 35/9, 10/9), 3, byrow=TRUE)
#' strategies <- c("Hawk", "Dove", "Mixed ESS")
#' ESset(A, strategies)
#' }


ESset <- function(A, strategies = c("1", "2", "3"), floats = TRUE) {
  if(!is.matrix(A) || !is.numeric(A)) {
    stop("A must be a numeric matrix.")
  }
  else if(nrow(A) != 3 || ncol(A) != 3) {
    stop("A must be of size 3x3.")
  }
  else if(length(strategies) != 3) {
    stop("Number of strategies must be 3.")
  }

  # calculate start point of ESset
  allESS <- ESS(A)
  stratESS <- ESS(A[-3, -3])
  ESset <- c()
  
  if(is.null(stratESS) || is.null(allESS)) {
    stop("Cannot calculate evolutionarily stable sets for models that do not have
      a proper ESS.")
  }
  
  if(!is.null(allESS)) {
    ESset <- rbind(ESset, allESS)
  }

  # calculate end point of ESset
  for(i in 1:nrow(stratESS)) {
    for(j in 0:2) {
      newStratESS <- matrix(append(stratESS, 0, j), 1)

      if(identical(newStratESS, allESS)) {
        break
      }

      ep <- c()

      for(i in 1:3) {
        ep <- rbind(ep, sum(A[i,] * newStratESS))
      }

      if(abs(max(ep) - min(ep)) < 1e-05) {
        if(!is.null(newStratESS)) {
          ESset <- rbind(ESset, newStratESS)
        }
      }
    }
  }

  if(is.null(ESset)) {
    stop("Cannot calculate evolutionarily stable sets for models that do not have
      a proper ESS.")
  }

  # change number representation
  if(!floats) {
    ESset <- MASS::fractions(ESset)
  }

  refSimp <- triangle()$coords

  # convert barycentric to cartesion coordinates
  setData <- geometry::bary2cart(refSimp, ESset)
  x <- y <- NULL
  setData <- data.frame(x = setData[, 1], y = setData[, 2])

  times <- seq(0, 100, by = 0.01)
  parameters <- c(A)

  states <- c()

  # generate sample states
  for(j in 1:3) {
    for(i in seq(0.1, 1, 0.1)) {
      if(j == 1){
        s1 <- i - 0.01
        s2 <- 0.01
        s3 <- 1 - s1 - s2
      }
      else if(j == 2) {
        s2 <- i - 0.01
        s3 <- 0.01
        s1 <- 1 - s2 - s3
      }
      else if(j == 3) {
        s3 <- i - 0.01
        s1 <- 0.01
        s2 <- 1 - s3 - s1
      }

      if(sum(s1, s2, s3) > 1 || sum(s1, s2, s3) < 0) {
        break
      }

      states <- rbind(states, c(s1, s2, s3))
    }
  }

  # create trajectories out of sample states
  odeData <- c()

  for(i in 1:nrow(states)) {
    out <- deSolve::ode(y = states[i, ], times = times, func = Replicator,
                        parms = parameters)
    odeData <- rbind(odeData, out[, -1])
  }

  odeData <- geometry::bary2cart(refSimp, odeData)
  odeData <- data.frame(x = odeData[, 1], y = odeData[, 2])

  # determine direction of trajectories
  step <- c()
  fac <- length(times) / nrow(states)

  for(i in seq(1, nrow(odeData) - fac, length(times))) {
    step <- c(step, i + fac)
  }

  step2 <- step + 1

  xend <- yend <- NULL
  arrowData <- data.frame(
    x = c(odeData[step, 1]),
    y = c(odeData[step, 2]),
    xend = c(odeData[step2, 1]),
    yend = c(odeData[step2, 2])
  )

  # prepare ESset plot
  trian <- triangle(strategies)$canvas +
    ggplot2::geom_line(
      data = setData,
      ggplot2::aes(x = x, y = y),
      size = 0.5
    ) +
    ggplot2::geom_point(
      data = odeData,
      ggplot2::aes(x = x, y = y),
      size = 0.1,
      shape = 16
    ) +
    ggplot2::geom_segment(
      data = arrowData,
      size = 0,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(0.25, "cm"), type = "closed"
      )
    )

  # plot ESset
  print(trian)

  return(ESset)
}
