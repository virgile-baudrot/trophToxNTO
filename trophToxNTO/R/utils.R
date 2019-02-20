#' @title ggplot2 dual axes
#'
#' @description ggplot2 dual axis
#'
#' @details see https://gist.github.com/jslefche/e4c0e9f57f0af49fca87
#'
#' @importFrom gtable gtable_add_grob
#' @importFrom grid grid.newpage
#' @import ggplot2
#'
#' @examples
#' # Create fake data.frame
#' data.add.x = data.frame(
#'   y1 = runif(100, 0, 100),
#'   x1 = runif(100, 0, 100)
#' )
#' # Add second x-axis that scales with first
#' data.add.x$x2 = (data.add.x$x1 + 50)^0.75
#'
#' # Create plots
#' plot1.x = qplot(y = y1, x = x1, data = data.add.x)
#' plot2.x = qplot(y = y1, x = x2, data = data.add.x)
#'
#' # Run function
#' ggplot_dual_axis(plot1.x, plot2.x, "x")
#'
#' @export
#'
ggplot_dual_axis = function(plot1, plot2, which.axis = "x") {

  plot2 <- plot2 + theme(axis.title.y = element_text(angle = 270))
  # Update plot with transparent panel
  #plot2 = plot2 + theme(panel.background = element_rect(fill = NA))

  grid.newpage()

  # Increase right margin if which.axis == "y"
  if(which.axis == "y") plot1 = plot1 + theme(plot.margin = unit(c(0.7, 1.5, 0.4, 0.4), "cm"))

  # Extract gtable
  g1 = ggplot_gtable(ggplot_build(plot1))

  g2 = ggplot_gtable(ggplot_build(plot2))

  # Overlap the panel of the second plot on that of the first
  pp = c(subset(g1$layout, name == "panel", se = t:r))

  g = gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)

  # Steal axis from second plot and modify
  axis.lab = ifelse(which.axis == "x", "axis-b", "axis-l")

  ia = which(g2$layout$name == axis.lab)

  ga = g2$grobs[[ia]]

  ax = ga$children[[2]]

  # Switch position of ticks and labels
  if(which.axis == "x") ax$heights = rev(ax$heights) else ax$widths = rev(ax$widths)

  ax$grobs = rev(ax$grobs)

  if(which.axis == "x")

    ax$grobs[[2]]$y = ax$grobs[[2]]$y - unit(1, "npc") + unit(0.15, "cm") else

      ax$grobs[[1]]$x = ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")

  # Modify existing row to be tall enough for axis
  if(which.axis == "x") g$heights[[2]] = g$heights[g2$layout[ia,]$t]

  # Add new row or column for axis label
  if(which.axis == "x") {

    g = gtable_add_grob(g, ax, 2, 4, 2, 4)

    g = gtable_add_rows(g, g2$heights[1], 1)

    g = gtable_add_grob(g, g2$grob[[6]], 2, 4, 2, 4)

  } else {

    g = gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)

    g = gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

    # g = gtable_add_grob(g, g2$grob[[7]], pp$t, length(g$widths), pp$b - 1)
    g = gtable_add_grob(g, g2$grob[[13]], pp$t, length(g$widths), pp$b - 1)

  }

  # Draw it
  grid.draw(g)

}



#' tag_facet
#'
#' @description Adds a dummy text layer to a ggplot to label facets and sets facet strips to blank.
#' This is the typical formatting for some journals that consider facets as subfigures
#' and want to minimise margins around figures.
#'
#' @details: see: https://github.com/baptiste/egg
#'
#' @param p ggplot
#' @param open opening character, default: (
#' @param close  closing character, default: )
#' @param tag_pool character vector to pick tags from
#' @param x x position within panel, default: -Inf
#' @param y y position within panel, default: Inf
#' @param hjust hjust
#' @param vjust vjust
#' @param fontface fontface
#' @param family font family
#' @param ... further arguments passed to geom_text layer
#'
#' @return plot with facet strips removed and replaced by in-panel tags
#' @importFrom ggplot2 geom_text ggplot_build theme element_blank aes
#' @export
#'
#' @examples
#' library(ggplot2)
#' mydf = data.frame(
#'   x = 1:90,
#'   y = rnorm(90),
#'   red = rep(letters[1:3], 30),
#'   blue = c(rep(1, 30), rep(2, 30), rep(3, 30)))
#'
#' p <- ggplot(mydf) +
#'   geom_point(aes(x = x, y = y)) +
#'   facet_wrap(
#'     ~ red + blue)
#' tag_facet(p)
#'
tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf,
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {

  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)

  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust,
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) + theme(strip.text = element_blank(),
                                                                                                  strip.background = element_blank())
}
