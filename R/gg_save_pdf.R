
#' Title
#'
#' @param list blabla
#' @param path blabla
#' @param filename blabla
#'
#' @return blabla
#'
#' @examples
#' blabla
#'
#'
gg_save_pdf = function(list, path, filename) {
  setwd(path)
  pdf(filename)
  for (p in list) {
    print(p)
  }
  dev.off()
  invisible(NULL)
}


