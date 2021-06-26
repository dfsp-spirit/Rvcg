
#' @title Compute one-ring vertex neigbors of vertices.
#'
#' @param x tmesh3d instance
#'
#' @param vi optional, vector of positive vertex indices for which to compute neighborhood, all are used if left at \code{NULL}.
#'
#' @return list of integer vectors, the neighborhoods.
#'
#' @examples
#' data(humface)
#' neighborhoods <- vcgOneRing(humface)
#' \dontrun{
#'   if(requireNamespace("fsbrain", quitely=TRUE)) {
#'   sjd = fsaverage.path(TRUE);
#'   surface = subject.surface(sjd, 'fsaverage', surface = "white", hemi = "lh");
#'   neigh = vcgOneRing(fs.surface.to.tmesh3d(surface));
#'   fsbrain::highlight.vertices.on.subject(sjd, 'fsaverage', verts_lh=neigh[[100]]);
#'   # assert max(unlist(lapply(neigh, max))) == 1643842 }
#'
#' @export
vcgOneRing <- function(x, vi=NULL, numstep=1L) {
  if(is.null(vi)) {
    vi = seq(ncol(x$vb));
  }
  vi <- as.integer(vi - 1L)
  vb <- x$vb
  it <- x$it - 1L
  out <- .Call("RVVadj",vb,it,vi,numstep)
}
