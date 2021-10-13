

#' @title Smooth per-vertex data on a mesh.
#'
#' @param x tmesh instance, the source mesh.
#'
#' @param data numeric vector, one value per mesh vertex.
#'
#' @param fwhm scalar double smoothing kernel full width at half max
#'
#' @param trunc_factor scalar double, truncation factor for Gaussian neighborhood, in Gaussian standard deviations.
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

  if(! is.vector(data)) {
    stop("Parameter 'data' must be a numeric vector.");
  }

  if(length(data) != ncol(x$vb)) {
    stop("Length of 'data' must match mesh 'x' vertex count.");
  }

  if(! is.numeric(fwhm)) {
    stop("Parameter 'fwhm' must be a scalar numeric value.");
  }
  if(! is.numeric(trunc_factor)) {
    stop("Parameter 'trunc_factor' must be a scalar numeric value.");
  }
  if(trunc_factor <= 0.0) {
    stop("Parameter 'trunc_factor' must be positive.");
  }
  if(fwhm <= 0.0) {
    stop("Parameter 'fwhm' must be positive.");
  }

  vb <- x$vb;
  it <- x$it - 1L;
  out <- .Call("RsmoothPVD",vb,it,data,fwhm,trunc_factor);
  return(out);
}



#' @title Internal testing function, do not use.
#'
#' @note This is NOT part of the API and may disappear/change at any time without notice, do NOT call this in your client code.
#'
#' @keywords internal
internal_test_with_brainmesh <- function() {
  warning("internal_test_with_brainmesh: This is NOT part of the API and may disappear/change at any time without notice, do NOT call this in your client code.");
  sj = "fsaverage5";
  #fsbrain::download_fsaverage3(T);
  sjd = fsbrain::fsaverage.path(T);
  lh = fsbrain::subject.surface(sjd, sj, hemi="lh");
  pvd = fsbrain::subject.morph.native(sjd, sj, "thickness", hemi="lh");
  fsbrain::vis.data.on.subject(sjd, sj, morph_data_lh = pvd);
  smoothed_pvd = Rvcg::vcgSmoothPVD(fsbrain::fs.surface.to.tmesh3d(lh), pvd, fwhm=5.0);
  fsbrain::vis.data.on.subject(sjd, sj, morph_data_lh = smoothed_pvd);
}


