#' @noRd
# save figure
save_figs <- function(name,
                      fig,
                      width = 6,
                      height = 6,
                      root = file.path(here::here(), "analysis/plots")) {

  dir.create(root, showWarnings = FALSE)
  fig_path <- function(name) {paste0(root, "/", name)}

  cowplot::save_plot(filename = fig_path(paste0(name,".png")),
                     plot = fig,
                     base_height = height,
                     bg = "white",
                     base_width = width)

  pdf(file = fig_path(paste0(name,".pdf")), width = width, height = height)
  print(fig)
  dev.off()

}
