

#' @title Smooth per vertex data on a mesh.
#' @param mesh triangular mesh of class \code{mesh3d}
#' @param data numeric vector of per vertex data, one scalar value per vertex
#' @param num_iter positive integer, the number of iterations
#' @param k the k for computing the k-ring-neighborhood. I.e., the max distance (in number of mesh edges) a vertex is allowed to be from the source vertex to be included in its neighborhood.
#' @return numeric vector, the smoothed data
#' @note If you want to smooth several data vectors on the same mesh, it is a lot faster to compute the neighborhood only once using \code{vcgVertexNeighbors} and then use \code{pvd_smoothnn_neigh} many times.
#' @examples
#'   data(humface);
#'   mean_curv = vcgCurve(humface)$meanvb;
#'   sm_mean_curv = pvd_smoothnn(humface, mean_curv, num_iter=50, k=3);
#'   hist(mean_curv); hist(sm_mean_curv);
#' @export
pvd_smoothnn <- function(mesh, data, num_iter, k=1L) {
  if(! (is.numeric(data) && is.vector(data))) {
    stop("Parameter 'data' must be a numeric vector.");
  }
  num_iter = as.integer(num_iter);
  if(length(num_iter) != 1L || num_iter < 1L) {
    stop("Parameter 'num_iter' must be a scalar, positive integer.");
  }
  k = as.integer(k);
  if(length(k) != 1L || k < 1L) {
    stop("Parameter 'k' must be a scalar, positive integer.");
  }

  neighborhood = Rvcg::vcgVertexNeighbors(mesh, numstep=k, include_self=TRUE);
  return(pvd_smoothnn_neigh(mesh, data, num_iter, neighborhood));
}


#' @title Perform NN smoothing in given neighborhood.
pvd_smoothnn_neigh <- function(mesh, data, num_iter, neighborhood) {
  if(! (is.numeric(data) && is.vector(data))) {
    stop("Parameter 'data' must be a numeric vector.");
  }

  num_iter = as.integer(num_iter);
  if(length(num_iter) != 1L || num_iter < 1L) {
    stop("Parameter 'num_iter' must be a scalar, positive integer.");
  }

  k = as.integer(k);
  if(length(k) != 1L || k < 1L) {
    stop("Parameter 'k' must be a scalar, positive integer.");
  }

  vb <- mesh$vb;
  it <- mesh$it - 1L;
  smoothed <- .Call("Rpvd_smoothnn", vb, it, data, num_iter, neighborhood);
  return(out);
}

