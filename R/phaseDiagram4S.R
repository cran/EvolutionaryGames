#' @name phaseDiagram4S
#' @title Phase Diagram for two-player games with four strategies
#' @description Plots phase diagram of a game with two players and four
#'  strategies.
#' @aliases phaseDiagram4S
#' @export phaseDiagram4S
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param A Numeric matrix of size 4x4 representing the number of strategies of
#'  a symmetric matrix game.
#' @param dynamic Function representing an evolutionary dynamic.
#' @param params Numeric vector with additional parameters for the evolutionary
#'  dynamic.
#' @param trajectory Numeric vector of size 4 representing the initial value
#'  for the trajectory to be examined.
#' @param strategies String vector of length 4 that names all strategies.
#' @param noRGL Logical value that handles diagram rotation. If
#'  set to \code{FALSE}, diagram will be rotatable, otherwise not. Default
#'  is \code{TRUE}.
#' @return None.
#' @examples
#' A <- matrix(c(5, -9, 6, 8, 20, 1, 2, -18, -14, 0, 2, 20, 13, 0, 4, -13),
#'  4, 4, byrow=TRUE)
#' state <- c(0.3, 0.2, 0.1, 0.4)
#' phaseDiagram4S(A, Replicator, NULL, state)

phaseDiagram4S <- function(A, dynamic, params = NULL, trajectory = NULL,
  strategies = c("1", "2", "3", "4"), noRGL = TRUE) {
  if(!is.matrix(A) || !is.numeric(A)) {
    stop("A must be a numeric matrix.")
  }
  else if(nrow(A) != 4 || ncol(A) != 4) {
    stop("A must be of size 4x4.")
  }
  else if(!is.null(params) && !is.numeric(params)){
    stop("params must be numeric.")
  }
  else if(!is.numeric(trajectory)){
    stop("trajectory must be a numeric vector.")
  }
  else if(length(trajectory) != 4){
    stop("trajectory must be of length 4.")
  }
  else if(!is.null(strategies) && length(strategies) != 4) {
    stop("Number of strategies does not match the number of columns of A.")
  }

  # group all parameters
  parameters <- as.vector(A)

  if(!is.null(params)) {
    p <- as.vector(params)
    parameters <- c(parameters, p)
  }

  times <- seq(0, 100, by = 0.01)

  # create trajectory
  out <- deSolve::ode(y = trajectory, times = times, func = dynamic,
                      parms = parameters)
  odeData <- out[, -1]
  
  # create rotatable diagram
  if(!noRGL) {
    x <- c(0, 0, 1, 1)
    y <- c(0, 1, 1, 0)
    z <- c(0, 1, 0, 1)
    
    for(i in 1:length(x)) {
      n <- setdiff(1:length(x), i)
      rgl::triangles3d(x[n], y[n], z[n], alpha = 0.5, col = "white")
      rgl::rgl.texts(x[i], y[i], z[i], text = strategies[i], col = "black")
    }   
    
    refSimp2 <- cbind(x, y, z)
    odeData <- geometry::bary2cart(refSimp2, odeData)
    rgl::rgl.points(odeData[,1], odeData[,2], odeData[,3], col = "black")
  }
  else {
    refSimp <- cbind(c(0.5, 0.5, 0, 1), c(sqrt(3)/2, sqrt(3)/6, 0, 0))
    # convert barycentric to cartesion coordinates
    odeData <- geometry::bary2cart(refSimp, odeData)
    odeData <- data.frame(x = odeData[, 1], y = odeData[, 2])

    ids <- factor(c("1", "2", "3"))
    id <- x <- y <- NULL
    
    # build 3-simplex
    coords <- data.frame(
      id = rep(ids, each = 3),
      x = c(refSimp[1, 1],
            refSimp[2, 1],
            refSimp[3, 1],
            refSimp[2, 1],
            refSimp[3, 1],
            refSimp[4, 1],
            refSimp[1, 1],
            refSimp[2, 1],
            refSimp[4, 1]),

      y = c(refSimp[1, 2],
            refSimp[2, 2],
            refSimp[3, 2],
            refSimp[2, 2],
            refSimp[3, 2],
            refSimp[4, 2],
            refSimp[1, 2],
            refSimp[2, 2],
            refSimp[4, 2])
    )

    value <- NULL
    vals <- data.frame(
      id = ids,
      value = c(0, 1, 2)
    )

    polyhData <- merge(coords, vals, by = c("id"))

    # prepare edge point labels
    labelsData <- data.frame(
      x = c(refSimp[1, 1], refSimp[2, 1],  refSimp[3, 1], refSimp[4, 1]),
      y = c(refSimp[1, 2] + 0.04,
            refSimp[2, 2] - 0.04,
            refSimp[3, 2] - 0.04,
            refSimp[4, 2] - 0.04)
    )

    # plot phase diagram
    ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = polyhData,
        ggplot2::aes(x = x, y = y, fill = value, group = id),
        color="black",
        size = 0.3
      ) +
      ggplot2::geom_point(
        data = odeData,
        ggplot2::aes(x = x, y = y),
        size = 0.1, shape = 16
      ) +
      ggplot2::scale_fill_gradient(low = "white", high = "gray93") +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::scale_x_continuous(limits = c(-0.2, 1.2)) +
      ggplot2::guides(fill = FALSE) +
      ggplot2::geom_text(
        data = labelsData,
        ggplot2::aes(x = x, y = y, label = strategies),
        size = 3.5
      )
  }
}
