
#' @title Compute one-ring vertex neigbors of vertices.
#'
#' @param x tmesh3d instance
#'
#' @param vi optional, vector of positive vertex indices for which to compute neighborhood, all are used if left at NULL.
#'
#' @return list of integer vectors, the neighborhoods.
#'
#' @export
vcgOneRing <- function(x, vi=NULL) {
  if(is.null(vi)) {
    vi = seq(ncol(x$vb));
  }
  vi <- as.integer(vi - 1L)
  vb <- x$vb
  it <- x$it - 1L
  out <- .Call("RVVadj",vb,it,vi)
  return(out)
}
