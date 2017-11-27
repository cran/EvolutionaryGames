#' @name phaseDiagram2S
#' @title Phase Diagram for two-player games with two strategies
#' @description Plots phase diagram of a game with two players and two
#'  strategies.
#' @aliases phaseDiagram2S
#' @export phaseDiagram2S
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param A Numeric matrix of size 2x2 representing the number of strategies of
#'  a symmetric matrix game.
#' @param dynamic Function representing an evolutionary dynamic.
#' @param params Numeric vector representing additional parameters for the
#'  evolutionary dynamic.
#' @param vectorField Logical value that handles vector field presentation. If
#'  set to \code{TRUE}, vector field will be shown, otherwise not. Default is
#'  \code{TRUE}.
#' @param strategies String vector of length 2 that names all strategies.
#' @return None.
#' @examples
#' A <- matrix(c(-1, 4, 0, 2), 2, 2, byrow=TRUE)
#' phaseDiagram2S(A, Replicator, strategies = c("Hawk", "Dove"))

phaseDiagram2S <- function(A, dynamic, params = NULL, vectorField = TRUE, strategies = c("1", "2")) {
  if(!is.matrix(A) || !is.numeric(A)) {
    stop("A must be a numeric matrix.")
  }
  else if(nrow(A) != 2 || ncol(A) != 2) {
    stop("A must be of size 2x2.")
  }
  else if(!is.null(params) && !is.numeric(params)){
    stop("params must be numeric.")
  }

  # group all parameters
  param <- as.vector(A)

  if(!is.null(params)) {
    param <- as.vector(params)
    param <- c(param, param)
  }

  times <- seq(0, 1, 0.01)
  y <- c()

  # obtain sample values
  for(i in times) {
    y <- c(y, dynamic(state = c(i, 1 - i), parameters = param)[[1]][1])
  }

  dynData <- data.frame(x = times, y = y)

  # create vector field
  num <- 10
  dist <- 1 / (num + 1)

  step <- seq(dist, 1 - dist, dist)
  x <- xend <- step

  for(i in 1:num) {
    s <- step[i]

    val <- dynamic(state = c(s, 1 - s), parameters = param)[[1]][1]
    fac <- 0.001

    if(val > 0) {
      xend[i] <- xend[i] + fac
    }
    else if(val < 0) {
      x[i] <- x[i] + fac
    }
  }

  arrData <- data.frame(x = x, xend = xend)
  lineData <- data.frame(x = c(0, 1), y = c(0, 0))
  xAxis <- paste("population share of strategy", strategies[1])

  # prepare phase diagram
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = dynData,
      ggplot2::aes(x = times, y = y),
      size = 0.3,
      color = "black"
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = xAxis, y = "dx/dt") +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, 1, 1),"cm"))

  vField <- list(
    ggplot2::geom_segment(
      data = arrData,
      size = 0.3,
      ggplot2::aes(x = x, xend = xend, y = 0,  yend = 0),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"))),
    ggplot2::geom_line(data = lineData, ggplot2::aes(x = x, y = y), size = 0.3)
  )

  # plot phase diagram
  if(vectorField) {
    print(p + vField)
  }
  else {
    print(p)
  }
}
