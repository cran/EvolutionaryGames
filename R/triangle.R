#' @name triangle
#' @title Triangle for 2-simplex operations
#' @description Generates a triangle representing the 2-simplex.
#' @aliases triangle
#' @export triangle
#' @author Daniel Gebele \email{dngebele@@gmail.com}
#' @param labels String vector of length 3 that names the edges of the
#'  triangle.
#' @return List of size 2 with members \code{coords} and \code{canvas}.
#'  \code{coords} holds edge coordinates of the 2-simplex, \code{canvas}
#'  a ggplot2 plot object of the 2-simplex.
#' @examples
#' triangle()

triangle <- function(labels = c("1", "2", "3")) {
  if(length(labels) != 3) {
    stop("Number of lables must be 3.")
  }

  # edge coordinates
  refSimp <- cbind(c(0.5, 0, 1), c(sqrt(3)/2, 0, 0))

  # prepare edge point labels
  x <- y <- NULL
  triangleData <- data.frame(x = refSimp[, 1], y = refSimp[, 2])
  labelsData <- data.frame(
    x = refSimp[, 1],
    y = refSimp[, 2] + c(0.04, -0.04, -0.04)
  )

  # prepare triangle plot
  plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = triangleData,
      ggplot2::aes(x = x, y = y),
      color = "black",
      fill = "white",
      size = 0.3
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::scale_x_continuous(limits = c(-0.2, 1.2)) +
    ggplot2::geom_text(
      data = labelsData,
      ggplot2::aes(x = x, y = y, label = labels),
      size = 3.5
    )

  return(list(coords = refSimp, canvas = plot))
}
