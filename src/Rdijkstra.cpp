#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include<vcg/complex/algorithms/geodesic.h>
#include <cassert>
#include <iostream>

using namespace tri;
using namespace Rcpp;

RcppExport SEXP Rdijkstra(SEXP vb_, SEXP it_, SEXP verts_, SEXP maxdist_)
{
  try {
    // Declare Mesh and helper variables
    IntegerVector verts(verts_);
    int n = verts.length();
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    double maxdist = as<double>(maxdist_);

    // Allocate mesh and fill it
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    m.vert.EnableVFAdjacency();
    m.vert.EnableQuality();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    tri::UpdateTopology<MyMesh>::VertexFace(m);

    // Prepare seed vector
    std::vector<MyVertex*> seedVec;
    for (int i=0; i < n; i++) {
      vi = m.vert.begin()+verts[i];
      seedVec.push_back(&*vi);
    }

    std::vector<MyVertex*> *inInterval;
    typename MyMesh::template PerVertexAttributeHandle<VertexPointer> sourcesHandle;
    sourcesHandle =  tri::Allocator<MyMesh>::AddPerVertexAttribute<MyMesh::VertexPointer> (m,"sources");
    typename MyMesh::template PerVertexAttributeHandle<VertexPointer> parentHandle;
    parentHandle =  tri::Allocator<MyMesh>::AddPerVertexAttribute<MyMesh::VertexPointer> (m,"parent");


    // Compute pseudo-geodesic distance by summing dists along shortest path in graph.
    tri::EuclideanDistance<MyMesh> ed;
    tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed, maxdist, inInterval, &sourcesHandle, &parentHandle);
    std::vector<float> geodist;
    vi=m.vert.begin();
    for (int i=0; i < m.vn; i++) {
      geodist.push_back(vi->Q());
      ++vi;
    }
    //return List::create(geodist); // also works
    return wrap(geodist);

  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RdijkstraPath(SEXP vb_, SEXP it_, SEXP vertsource_, SEXP verttarget_, SEXP maxdist_)
{
  try {
    // Declare Mesh and helper variables
    IntegerVector verts_source(vertsource_); // contains only 1 vertex index
    IntegerVector verts_target(verttarget_); // contains only 1 vertex index
    int n = verts_source.length();
    assert(n==1);
    assert(verts_target.length()==1);
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    double maxdist = as<double>(maxdist_);

    // Allocate mesh and fill it
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    m.vert.EnableVFAdjacency();
    m.vert.EnableQuality();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    tri::UpdateTopology<MyMesh>::VertexFace(m);

    // Prepare seed vector with a single vertex
    std::vector<MyVertex*> seedVec;
    vi = m.vert.begin()+verts_source[0];
    seedVec.push_back(&*vi);

    std::vector<MyVertex*> *inInterval;
    typename MyMesh::template PerVertexAttributeHandle<VertexPointer> sourcesHandle;
    sourcesHandle =  tri::Allocator<MyMesh>::AddPerVertexAttribute<MyMesh::VertexPointer> (m,"sources");
    typename MyMesh::template PerVertexAttributeHandle<VertexPointer> parentHandle;
    parentHandle =  tri::Allocator<MyMesh>::AddPerVertexAttribute<MyMesh::VertexPointer> (m,"parent");


    // Compute pseudo-geodesic distance by summing dists along shortest path in graph.
    tri::EuclideanDistance<MyMesh> ed;
    tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed, maxdist, inInterval, &sourcesHandle, &parentHandle);

    std::vector<float> geodist;
    vi=m.vert.begin();
    for (int i=0; i < m.vn; i++) {
      geodist.push_back(vi->Q());
      ++vi;
    }

    std::vector<int> geopath;
    // TODO: compute geodesic path

    return List::create(geodist, geopath);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

// Compute for each vertex the geodesic neighborhood with radius x, i.e., all vertices which are in geodesic distance <= x.
RcppExport SEXP Rgeodesicneigh(SEXP vb_, SEXP it_, SEXP neighdist_, SEXP ignore_mask_)
{
  try {
    // Declare Mesh and helper variables
    MyMesh m;
    VertexIterator vi, vo;
    FaceIterator fi;
    float neighdist = Rcpp::as<float>(neighdist_);

    IntegerVector ignore_mask(ignore_mask_);

    // Allocate mesh and fill it
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    m.vert.EnableVFAdjacency();
    m.vert.EnableQuality();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    tri::UpdateTopology<MyMesh>::VertexFace(m);

    std::vector<std::vector<int>> all_neighborhoods;
    for (int cur_vert=0; cur_vert < m.vn; cur_vert++) {
      // Prepare seed vector with single vertex.
      std::vector<MyVertex*> seedVec;
      vi = m.vert.begin()+cur_vert;
      seedVec.push_back(&*vi);


      // Compute pseudo-geodesic distance by summing dists along shortest path in graph.
      std::vector<int> cur_neighborhood;

      if(ignore_mask[cur_vert] == 0) { // Only compute for non-ignored vertices. This returns an empty neighborhood for ignored ones.
        tri::EuclideanDistance<MyMesh> ed;
        tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed);
        vo=m.vert.begin();
        for (int i=0; i < m.vn; i++) {
          if(vo->Q() < neighdist && ignore_mask[i] == 0) {  // Only add non-ignored vertices to neighborhood.
            cur_neighborhood.push_back(i);
          }
          ++vo;
        }
      }
      all_neighborhoods.push_back(cur_neighborhood);
    }
    return wrap(all_neighborhoods);

  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

// Compute mean geodesic dist to all others for each vertex.
RcppExport SEXP Rgeodesicmeandist(SEXP vb_, SEXP it_, SEXP ignore_mask_)
{
  try {
    // Declare Mesh and helper variables
    MyMesh m;
    VertexIterator vi, vo;
    FaceIterator fi;
    IntegerVector ignore_mask(ignore_mask_);
    int nv = ignore_mask.length();

    int num_non_ignored_verts = 0;
    for (int i=0; i < nv; i++) {
      if(ignore_mask[i] == 0) {
        num_non_ignored_verts += 1;
      }
    }

    // Allocate mesh and fill it
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    m.vert.EnableVFAdjacency();
    m.vert.EnableQuality();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    tri::UpdateTopology<MyMesh>::VertexFace(m);

    std::vector<float> meangeodist;
    for (int cur_vert=0; cur_vert < m.vn; cur_vert++) {

      if(ignore_mask[cur_vert] == 1) { // do not compute mean distance for ignored vertices.
        meangeodist.push_back(-1);
      } else {
        // Compute pseudo-geodesic distance by summing Euclidean dists along shortest path in graph.
        tri::EuclideanDistance<MyMesh> ed;

        std::vector<MyVertex*> seedVec;
        vi = m.vert.begin()+cur_vert;
        seedVec.push_back(&*vi);

        tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed);
        double sumdist = 0.0;
        vo=m.vert.begin();
        for (int i=0; i < m.vn; i++) {
          if(ignore_mask[cur_vert] == 0) {
            sumdist += (vo->Q() / num_non_ignored_verts);
          }
          ++vo;
        }
        meangeodist.push_back(sumdist);
      }
    }
    return wrap(meangeodist);

  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

