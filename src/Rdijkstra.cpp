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

    int source_vertex_index = verts_source[0];
    int target_vertex_index = verts_target[0];
    std::vector<int> geopath;
    //(*vi).P().X()
    // TODO: compute geodesic path
    //geopath.push_back(*(sourcesHandle[]))
    std::cout << "Source vertex is " << source_vertex_index << ", target is " << target_vertex_index << "\n";
    //std::cout << "Length of sourcesHandle: " << sourcesHandle.size() << "\n";
    //std::cout << "Length of parentHandle: " << parentHandle.size() << "\n";
    //geopath.push_back(sourcesHandle[target_vertex_index]); // TODO: howto get vertex index from vertex pointer? here we push tje pointer (which is invalid to push to an int array)

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


RcppExport SEXP RGeodesicPath(SEXP vb_, SEXP it_, SEXP source_, SEXP targets_, SEXP maxdist_)
{
  try {
    // Declare Mesh and helper variables
    int source = Rcpp::as<int>(source_);
    IntegerVector targets(targets_);
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
    vi = m.vert.begin()+source;
    seedVec.push_back(&*vi);

    std::vector<MyVertex*> *inInterval;
    typename MyMesh::template PerVertexAttributeHandle<VertexPointer> sourcesHandle;
    sourcesHandle =  tri::Allocator<MyMesh>::AddPerVertexAttribute<MyMesh::VertexPointer> (m, "sources");
    typename MyMesh::template PerVertexAttributeHandle<VertexPointer> parentHandle;
    parentHandle =  tri::Allocator<MyMesh>::AddPerVertexAttribute<MyMesh::VertexPointer> (m, "parent");

    // Compute pseudo-geodesic distance by summing dists along shortest path in graph.
    tri::EuclideanDistance<MyMesh> ed;
    tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed, maxdist, inInterval, &sourcesHandle, &parentHandle);
    std::vector<float> geodist;
    vi=m.vert.begin();
    for (int i=0; i < m.vn; i++) {
      geodist.push_back(vi->Q());
      ++vi;
    }

    // Geodesic path computation (the backtracking part of Dijkstra's Algorithm)
    // !!! See also: https://stackoverflow.com/questions/10667867/argmin-for-vectordouble-in-c
    //     A comment suggests we can substract iterators from one another to get indices, e.g.,
    //    "Since you want an index and you are working with vectors, you can then substract the resulting iterator from vec.begin() to get such index."
    // This way we could use the existing pointers provided by VCGLIB.
    std::vector<std::vector<int>> paths;
    <std::vector<float> path_lengths;
    for(int i=0; i<targets.size(), ++i) {
      int target_vertex = targets[i];
      int current_vertex = target_vertex;
      std::vector<int> path;
      path.push_back(current_vertex);
      std::vector<int> visited;
      float path_length = 0.0f;
      while(current_vertex != source) {
        visited.push_back(current_vertex);

        // TODO: Compute vertex neighborhood in VCGLIB
        std::vector<int> neigh = mesh.vertex.neighbors(surface, source_vertices = current_vertex)$vertices;

        std::vector<int> neigh_unvisited; // Keep only unvisited neighbors and track their distance to source.
        std::vector<float> neigh_unvisited_dists;
        for(int neighidx=0; neighidx < neigh.size(), ++neighidx) {
          int neighvert = neigh[neighidx];
          if(std::find(neigh.begin(), neigh.end(), neighvert) == neigh.end()) {
            neigh_unvisited.push_back(neighvert);
            neigh_unvisited_dists.push_back(geodist[neighvert]);
          }
        }
        int closest_to_source = std::distance(neigh_unvisited.begin(), std::min_element(neigh_unvisited_dists.begin(), neigh_unvisited_dists.end()));
        path_length += std::min(neigh_unvisited_dists.begin(), neigh_unvisited_dists.end());
        path.push_back(closest_to_source);
        current_vertex = closest_to_source;
      }

      std::reverse(path.begin(), path.end());
      paths.push_back(path);
      path_lengths.push_back(path_length);
    }

    List L = List::create(Named("paths") = paths , _["path_lengths"] = path_lengths, _["geodist"] = geodist);
    return L;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

