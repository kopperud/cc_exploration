##
library(ggplot2)

geom_stepribbon <- function(
  mapping     = NULL,
  data        = NULL,
  stat        = "identity",
  position    = "identity",
  na.rm       = FALSE,
  show.legend = NA,
  inherit.aes = TRUE, ...) {
  
  layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomStepribbon,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(na.rm = na.rm, ... )
  )
}

GeomStepribbon <- ggplot2::ggproto(
  "GeomStepribbon", GeomRibbon,
  
  extra_params = c("na.rm"),
  
  draw_group = function(data, panel_scales, coord, direction = "vh", include_final = FALSE, na.rm = FALSE) {
    n <- nrow(data)
    data <- as.data.frame(data)[order(data$x), ]
    
    if (direction == "vh") {
      xs <- rep(1:n, each = 2)[-2*n]
      ys <- c(1, rep(2:n, each = 2))
    } else if (direction == "hv") {
      xs <- c(1, rep(2:n, each = 2))
      ys <- rep(1:n, each = 2)[-2*n]
    } else {
      abort("Parameter `direction` is invalid.")
    }
    if(!include_final){
      xs <- tail(xs, n = -1)
      ys <- tail(ys, n = -1)
    }
    x <- data$x[xs]
    ymin <- data$ymin[ys]
    ymax <- data$ymax[ys]
    data_attr <- data[xs, setdiff(names(data), c("x", "ymin", "ymax"))]
    data <- ggplot2:::new_data_frame(c(list(x = x, ymin = ymin, ymax = ymax), data_attr))
    
    GeomRibbon$draw_group(data, panel_scales, coord, na.rm = FALSE)
    
  }
  
)