
#' @title Compute one-ring vertex neigbors of vertices.
#'
#' @param x tmesh3d instance
#'
#' @return list of integer vectors, the neighborhoods. Curently this function still returns FACE indices.
#'
#' @examples
#' data(humface)
#' neighborhoods <- vcgOneRing(humface)
#'
#' @export
vcgOneRing <- function(x) {
  vb <- x$vb
  it <- x$it - 1L
  out <- .Call("RVFadj",vb,it)
  return(out)
  # neigh = out
  # return(unique(as.integer(surface$faces[neigh[[1000]],])))
}
