
#' @title Compute pseudo-geodesic distances on a triangular mesh
#' @param x triangular mesh of class \code{mesh3d}
#' @param vertpointer integer: references indices of vertices on the mesh
#' @return returns a vector of shortest distances for each of the vertices to one of the vertices referenced in \code{vertpointer}
#' @examples
#' ## Compute geodesic distance between all mesh vertices and the first vertex of a mesh
#' data(humface)
#' humface <- vcgIsolated(vcgClean(humface,sel=0:6,iterate=TRUE))
#' geo <- vcgDijkstra(humface,1)
#' if (interactive()) {
#' require(Morpho);require(rgl)
#' meshDist(humface,distvec = geo)
#' spheres3d(vert2points(humface)[1,],col=2)
#' }
#' @note Make sure to have a clean manifold mesh. Note that this computes the length of the pseudo-geodesic path (following the edges) between the two vertices.
#' @export
vcgDijkstra <- function(x, vertpointer) {
    vertpointer <- as.integer(vertpointer-1)
    vb <- x$vb
    it <- x$it-1
    out <- .Call("Rdijkstra",vb,it,vertpointer)
    return(out)
}


#' @title Compute, for each mesh vertex, all neighbors within a given geodesic distance.
#'
#' @inheritParams vcgDijkstra
#'
#' @param dist double, a single scalar defining the max geodesic distance for the neighborhood.
#'
#' @param ignore_mask logical vector of length \code{|V|}, the number of vertices in the mesh. Each position must indicate whether the vertex should be ignored. Keep in mind that ignoring vertices may lead to a neighborhood consisting of several isolated patches.
#'
#' @return list of integer vectors, the neighbors. The length of the outer list equals the number of vertices in the mesh.
#'
#' @examples
#' \dontrun{
#'   fsbrain::download_fsaverage3(TRUE);
#'   sjd = fsbrain::fsaverage.path();
#'   sj = "fsaverage3";
#'   sf = fsbrain::subject.surface(sjd, sj, "white", "lh");
#'   mask = fsbrain::subject.mask(sjd, sj, hemi = "lh", invert_mask = FALSE);
#'   tm = fsbrain::fs.surface.to.tmesh3d(sf);
#'   neigh = Rvcg::vcgGeodesicNeigh(tm, 15.0, ignore_mask = mask);
#'   fsbrain::highlight.vertices.on.subject(sjd, sj,
#'     verts_lh = neigh[[500]]); # show vertex 638 neighborhood
#' }
#'
#' @export
vcgGeodesicNeigh <- function(x, dist, ignore_mask = NULL) {
    if(is.null(ignore_mask)) {
      ignore_mask = rep(FALSE, dim(x$vb)[2]);
    }
    ignore_mask = as.integer(ignore_mask);
    num_verts = dim(x$vb)[2];
    if(length(ignore_mask) != num_verts) {
        stop(sprintf("Ignore mask length (%d) must equal vertex count in mesh (%d).\n", length(ignore_mask), num_verts));
    }
    if(sum(ignore_mask) == num_verts) {
        stop("Cannot ignore all vertices.");
    }
    vb <- x$vb;
    it <- x$it - 1L;
    out <- .Call("Rgeodesicneigh", vb, it, dist, ignore_mask);
    return(out);
}


#' @title Compute, for each mesh vertex, the mean geodesic distance to all other vertices.
#'
#' @inheritParams vcgGeodesicNeigh
#'
#' @return vector of doubles, the mean geodesic distances. The length of the vector equals the number of vertices in the mesh. Ignored vertices will have an \code{NA} value.
#'
#' @examples
#' \dontrun{
#'   fsbrain::download_fsaverage3(TRUE);
#'   sjd = fsbrain::fsaverage.path();
#'   sj = "fsaverage3";
#'   sf = subject.surface(sjd, sj, "white", "lh");
#'   mask = fsbrain::subject.mask(sjd, sj, hemi = "lh", invert_mask = FALSE);
#'   tm = fsbrain::fs.surface.to.tmesh3d(sf);
#'   md = Rvcg::vcgGeodesicMeanDist(tm, ignore_mask = mask);
#'   fsbrain::vis.data.on.subject(sjd, sj, morph_data_lh = md);
#' }
#'
#' @export
vcgGeodesicMeanDist <- function(x, ignore_mask = NULL) {
    if(is.null(ignore_mask)) {
      ignore_mask = rep(FALSE, dim(x$vb)[2]);
    }
    ignore_mask = as.integer(ignore_mask);
    num_verts = dim(x$vb)[2];
    if(length(ignore_mask) != num_verts) {
        stop(sprintf("Ignore mask length (%d) must equal vertex count in mesh (%d).\n", length(ignore_mask), num_verts));
    }
    if(sum(ignore_mask) == num_verts) {
        stop("Cannot ignore all vertices.");
    }
    vb <- x$vb;
    it <- x$it - 1L;
    out <- .Call("Rgeodesicmeandist", vb, it, ignore_mask);
    out[out < 0] = NA;
    return(out);
}


#' @title Compute pseudo-geodesic distance between two points on a mesh
#' @param x triangular mesh of class \code{mesh3d}
#' @param pt1 3D coordinate on mesh or index of vertex
#' @param pt2 3D coordinate on mesh or index of vertex
#' @return returns the geodesic distance between \code{pt1} and \code{pt2}.
#' @note Make sure to have a clean manifold mesh. Note that this computes the length of the pseudo-geodesic path (following the edges) between the two vertices closest to these points.
#' @examples
#' data(humface)
#' pt1 <- humface.lm[1,]
#' pt2 <- humface.lm[5,]
#' vcgGeodist(humface,pt1,pt2)
#' @export
vcgGeodist <- function(x,pt1,pt2) {
    if (length(pt1) == 1)
        pt1 <- x$vb[1:3,pt1]
    if (length(pt2) == 1)
        pt2 <- x$vb[1:3,pt2]
    mypts <- rbind(pt1,pt2)
    clost <- vcgKDtree(x,mypts,k=1)
    geo <- vcgDijkstra(x,vertpointer = clost$index[1,1])[clost$index[2,1]]
    return(geo)
}
