
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
#' @return list of integer vectors, the neighbors. The length of the outer list equals the number of vertices in the mesh.
#'
#' @examples
#' \dontrun{
#'   fsbrain::download_fsaverage3(TRUE);
#'   sjd = fsbrain::fsaverage.path();
#'   sf = subject.surface(sjd, "fsaverage3", "white", "lh");
#'   tm = fsbrain::fs.surface.to.tmesh3d(sf);
#'   neigh = Rvcg::vcgGeodesicNeigh(tm, 15.0);
#'   fsbrain::highlight.vertices.on.subject(sjd, "fsaverage3",
#'     verts_lh = neigh[[638]]); # show vertex 638 neighborhood
#' }
#'
#' @export
vcgGeodesicNeigh <- function(x, dist) {
    vb <- x$vb;
    it <- x$it - 1L;
    out <- .Call("Rgeodesicneigh", vb, it, dist);
    return(out);
}


#' @title Compute, for each mesh vertex, the mean geodesic distance to all other vertices.
#'
#' @inheritParams vcgDijkstra
#'
#' @return vector of doubles, the mean geodesic distances. The length of the vector equals the number of vertices in the mesh.
#'
#' @examples
#' \dontrun {
#'   fsbrain::download_fsaverage3(TRUE);
#'   sjd = fsbrain::fsaverage.path();
#'   sf = subject.surface(sjd, "fsaverage3", "white", "lh");
#'   tm = fsbrain::fs.surface.to.tmesh3d(sf);
#'   md = Rvcg::vcgGeodesicMeanDist(tm);
#'   fsbrain::vis.data.on.subject(sjd, "fsaverage3", morph_data_lh = md);
#' }
#'
#' @export
vcgGeodesicMeanDist <- function(x) {
    vb <- x$vb;
    it <- x$it - 1L;
    out <- .Call("Rgeodesicmeandist", vb, it);
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
