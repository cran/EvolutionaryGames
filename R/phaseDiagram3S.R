#' @name phaseDiagram3S
#' @title Phase Diagram for two-player games with three strategies
#' @description Plots phase diagram of a game with two players and three
#'  strategies.
#' @aliases phaseDiagram3S
#' @export phaseDiagram3S
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param A Numeric matrix of size 3x3 representing the number of strategies of
#'  a symmetric matrix game.
#' @param dynamic Function representing an evolutionary dynamic.
#' @param params Numeric vector with additional parameters for the evolutionary
#'  dynamic.
#' @param trajectories Numeric matrix of size mx3. Each row represents the
#'  initial values for the trajectory to be examined.
#' @param contour Logical value that handles contour diagram presentation. If
#'  set to \code{TRUE}, contour diagram will be shown, otherwise not. Default
#'  is \code{FALSE}.
#' @param vectorField Logical value that handles vector field presentation. If
#'  set to \code{TRUE}, vector field will be shown, otherwise not. Default is
#'  \code{FALSE}.
#' @param strategies String vector of length 3 that names all strategies.
#' @return None.
#' @examples
#' A <- matrix(c(0, -2, 1, 1, 0, -2, -2, 1, 0), 3, byrow=TRUE)
#' state <- matrix(c(0.4, 0.3, 0.3), 1, 3, byrow=TRUE)
#'
#' phaseDiagram3S(A, Replicator, NULL, state, FALSE, FALSE)
#' phaseDiagram3S(A, Replicator, NULL, state, TRUE, TRUE)
#' 
#' # Plot two trajectories rather than only one:
#' A <- matrix(c(0, -2, 1, 1, 0, -2, -2, 1, 0), 3, byrow=TRUE)
#' state <- matrix(c(0.4, 0.3, 0.3, 0.6, 0.2, 0.2), 2, 3, byrow=TRUE)
#' phaseDiagram3S(A, Replicator, NULL, state, FALSE, FALSE)

phaseDiagram3S <- function(A, dynamic, params = NULL, trajectories = NULL,
  contour = FALSE, vectorField = FALSE, strategies = c("1", "2", "3")) {
  if(!is.matrix(A) || !is.numeric(A)) {
    stop("A must be a numeric matrix.")
  }
  else if(nrow(A) != 3 || ncol(A) != 3) {
    stop("A must be of size 3x3.")
  }
  else if(!is.null(params) && !is.numeric(params)){
    stop("params must be numeric.")
  }
  else if(!is.null(trajectories) && !is.numeric(trajectories)){
    stop("trajectories must be a numeric matrix.")
  }
  else if(!is.null(trajectories) && ncol(trajectories) != 3){
    stop("trajectories must be of size mx3.")
  }
  else if(!is.null(strategies) && length(strategies) != 3) {
    stop("Number of strategies does not match the number of columns of A.")
  }
  
  # group all parameters
  parameters <- as.vector(A)
  x <- y <- NULL
  
  if(!is.null(params)) {
    p <- as.vector(params)
    parameters <- c(parameters, p)
  }
  
  times <- seq(0, 100, by = 0.01)
  refSimp <- triangle()$coords
  odeData <- arrowData <- c()
  
  # create trajectories
  if(!is.null(trajectories)) {
    for(i in 1:nrow(trajectories)) {
      out <- deSolve::ode(y = trajectories[i, ], times = times, func = dynamic,
                          parms = parameters)
      odeData <- rbind(odeData, out[, -1])
    }
    
    # convert barycentric to cartesian coordinates
    odeData <- geometry::bary2cart(refSimp, odeData)
    odeData <- data.frame(x = odeData[, 1], y = odeData[, 2])
    
    # determine direction of trajectories
    arrNum <- 20
    dist <- length(times) / (arrNum + 1)
    step <- c()
    
    for(i in seq(1, nrow(odeData), length(times))) {
      for(j in 1:arrNum) {
        step <- c(step, i + j * dist)
      }
    }
    
    step2 <- step + 1
    xend <- yend <- NULL
    arrowData <- data.frame(
      x = c(odeData[step, 1]),
      y = c(odeData[step, 2]),
      xend = c(odeData[step2, 1]),
      yend = c(odeData[step2, 2])
    )
  }
  
  maxVelocity <- 0
  density <- c()
  x <- y <- z <- c()
  
  # determine velocity
  for(i in seq(0, 1, by = 0.1)) {
    for(j in seq(0, 1, by = 0.1)) {
      if(i + j > 1) {
        break
      }
      
      x <- c(x, i)
      y <- c(y, j)
      z <- c(z, 1 - i - j)
      
      dX <- dynamic(state = c(i, j, 1 - i - j), parameters = parameters)[[1]]
      dist <- sqrt(dX[1]^2 + dX[2]^2 + dX[3]^2)
      
      density <- c(density, dist)
      maxVelocity <- max(dist, maxVelocity)
    }
  }
  
  contourData <- cbind(x, y, z)
  
  x1 <- y1 <- z1 <- c()
  x2 <- y2 <- z2 <- c()
  
  # calculate vector field
  for(i in seq(0, 1, by = 0.05)) {
    for(j in seq(0, 1, by = 0.05)) {
      if(i + j > 1) {
        break
      }
      
      x1 <- c(x1, i)
      y1 <- c(y1, j)
      z1 <- c(z1, 1 - i - j)
      
      dX <- dynamic(state = c(i, j, 1 - i - j), parameters = parameters)[[1]]
      dist <- sqrt(dX[1]^2 + dX[2]^2 + dX[3]^2)
      
      x2 <- c(x2, dX[1] * dist)
      y2 <- c(y2, dX[2] * dist)
      z2 <- c(z2, dX[3] * dist)
    }
  }
  
  vecData1 <- cbind(x1, y1, z1)
  vecData2 <- cbind(x2, y2, z2)
  density <- density / maxVelocity
  
  # convert barycentric to cartesion coordinates
  contourData <- geometry::bary2cart(refSimp, contourData)
  vecData1 <- geometry::bary2cart(refSimp, vecData1)
  vecData2 <- geometry::bary2cart(refSimp, vecData2)
  
  vecData <- cbind(vecData1, vecData2)
  vecData <- data.frame(
    x = vecData[, 1],
    y = vecData[, 2],
    xend = vecData[, 1] + vecData[, 3] * 0.5,
    yend = vecData[, 2] + vecData[, 4] * 0.5
  )
  
  if(any(vecData < 0)) {
    vecData$xend = vecData$x
    vecData$yend = vecData$y
  }
  
  resol <- 300
  
  # perform interpolation
  contourData <- interp::interp(
    contourData[, 1],
    contourData[, 2],
    density,
    seq(0, 1, length = resol),
    seq(0, 1, length = resol)
  )
  
  # convert data into long-format
  contourData <- reshape2::melt(contourData$z)
  
  contourData[, 1:2] <- contourData[, 1:2] / resol
  contourData <- data.frame(
    x = contourData[, 1],
    y = contourData[, 2],
    z = contourData[, 3]
  )
  
  # prepare phase diagram
  p <- triangle(strategies)$canvas
  pal <- grDevices::colorRampPalette(c("blue","cyan", "green", "yellow", "red"))(5)
  
  if(contour) {
    p <- p +
      ggplot2::geom_tile(
        data = contourData,
        ggplot2::aes(x = x, y = y, fill = z)
      ) +
      ggplot2::scale_fill_gradientn(colours = pal, na.value = NA) + 
      ggplot2::labs(fill = "Velocity")
  }
  
  if(!is.null(trajectories)) {
    p <- p +
      ggplot2::geom_point(
        data = odeData,
        ggplot2::aes(x = x, y = y),
        size = 0.1,
        shape = 16
      ) +
      ggplot2::geom_segment(
        data = arrowData,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        arrow = ggplot2::arrow(
          length = ggplot2::unit(0.25, "cm"),
          type = "closed"
        ),
        size = 0
      )
  }
  
  if(vectorField) {
    p <- p + ggplot2::geom_segment(
      data = vecData,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(0.1, "cm"),
        type = "closed"
      ),
      size = 0.4
    )
  }
  
  # plot phase diagram
  print(p)
}
