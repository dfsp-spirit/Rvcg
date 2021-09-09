

#' @title Smooth per-vertex data on mesh.
#'
#' @param x tmesh instance, the source mesh.
#'
#' @param data numeric vector, one value per vertex
#'
#' @param fwhm scalar double smoothing kernel full width at half max
#'
#' @param trunc_factor scalar double, truncation factor for Gaussian neighborhood
#'
#' @return numeric vector, the smoothed data
#'
#' @examples
#' data(humface)
#' pvd = rnorm(ncol(humface$vb), mean = 5.0, sd = 1.0);
#' smoothed_pvd = vcgSmoothPVD(humface, pvd, fwhm=5.0);
#'
#' @export
vcgSmoothPVD <- function(x, data, fwhm, trunc_factor=3.5) {

  vb <- x$vb;
  it <- x$it - 1L;
  out <- .Call("RsmoothPVD",vb,it,data,fwhm,trunc_factor);
  return(out);
}
