

#' @title Smooth per-vertex data on mesh.
#'
#' @param x tmesh instance, the source mesh.
#'
#' @param data numeric vector, one value per vertex
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

  vb <- x$vb;
  it <- x$it - 1L;
  out <- .Call("RsmoothPVD",vb,it,data,fwhm,trunc_factor);
  return(out);
}



#test_with_brainmesh <- function() {
#  fsbrain::download_fsaverage3();
#  sjd = fsbrain::fsaverage.path(T);
#  lh = fsbrain::subject.surface(sjd, "fsaverage3", hemi="lh");
#  pvd = fsbrain::subject.morph.native(sjd, "fsaverage3", "thickness", hemi="lh");
#  fsbrain::vis.data.on.subject(sjd, "fsaverage3", morph_data_lh = pvd);
#  smoothed_pvd = vcgSmoothPVD(fsbrain::fs.surface.to.tmesh3d(lh), pvd, fwhm=5.0);
#  fsbrain::vis.data.on.subject(sjd, "fsaverage3", morph_data_lh = smoothed_pvd);
#}

