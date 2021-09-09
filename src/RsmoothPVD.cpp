
#include "typedef.h"
#include "RvcgIO.h"

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <cassert>
//#include <iostream>


std::vector<float> dijkstra(MyMesh& m, std::vector<int> verts, float maxdist); // declaration for function from Rdijkstra.cpp

/// Used in C++ code only, not exported.
inline float fhwm_to_gstd(const float fwhm) {
  return(fwhm / sqrt(log(256.0)));
}


/// Compute Gaussian weights for the neighborhoods of all vertices, based on geodesic distances.
/// Used in C++ code only, not exported.
std::vector<std::vector<float> > gauss_weights(const std::vector<std::vector<int> > geod_neigh_indices, const std::vector<std::vector<float> > geod_neigh_dists, const float gstd) {
  std::vector<std::vector<float> > weights(geod_neigh_indices.size());

  assert(geod_neigh_indices.size() == geod_neigh_dists.size());

  float gvar2 = 2 * (gstd * gstd);
  float f = 1.0 / (sqrt(2 * M_PI) * gstd);
  float gsum;
  for(size_t i=0; i<geod_neigh_indices.size(); i++) { // iterate over vertex count in mesh
    gsum = 0.0;
    std::vector<float> vertex_weights(geod_neigh_indices[i].size());
    size_t local_idx = 0L;
    for(size_t j=0; j<geod_neigh_indices[i].size(); j++) {
      float d = geod_neigh_dists[i][j];
      float g = f * exp(-(d * d) / (gvar2));
      vertex_weights[j] = g;
      gsum += g;
    }
    for(size_t j=0; j<geod_neigh_indices[i].size(); j++) {
      vertex_weights[j] /= gsum;
    }
    weights[i] = vertex_weights;
  }
  return(weights);
}


// Apply Gaussian weights to neighborhood data values to smooth data.
/// Used in C++ code only, not exported.
std::vector<float> spatial_filter(const std::vector<float> data, const std::vector<std::vector<int> > geod_neigh_indices, const std::vector<std::vector<float> > geod_neigh_gauss_weights) {
  std::vector<float> smoothed_data(data.size());
  float smoothed_val;
  for(size_t i=0; i<data.size(); i++) {
    smoothed_val = 0.0;
    for(size_t j=0; j<geod_neigh_indices[i].size(); j++) {
      smoothed_val += data[geod_neigh_indices[i][j]] * geod_neigh_gauss_weights[i][j];
    }
    smoothed_data[i] = smoothed_val;
  }
  return(smoothed_data);
}


/// Perform Gaussian smoothing of the given per-vertex data for the mesh.
/// @param vb_ the xyz vertex coordinates of the mesh
/// @param it_ the faces of the mesh, given as indices into the vertex list
/// @param data_ an R numerical vector with one value per mesh vertex
/// @param fwhm_ the FWHM for the Gaussian kernel
/// @param truncfactor_ the cutoff factor after which to end the Gaussian neighborhood, in Gaussian standard deviations
/// NOTE: This currently computes the full mesh neighborhood distances at once, which may result in out-of-memory issues for large meshes.
/// NOTE2: This function currently does not support NA values in the data.
RcppExport SEXP RsmoothPVD(SEXP vb_, SEXP it_, SEXP data_, SEXP fwhm_, SEXP truncfactor_) {
  float fwhm = Rcpp::as<float>(fwhm_);
  float gstd = fhwm_to_gstd(fwhm);
  float maxdist = gstd * Rcpp::as<float>(truncfactor_);
  std::vector<float> data = Rcpp::as<std::vector<float>>(data_);

  MyMesh m;
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);

  assert(m.vn == data.size());

  std::vector<std::vector<int> > geod_indices;
  std::vector<std::vector<float> > geod_distances;
  std::vector<int> verts(1);

  for(int i=0; i<data.size(); i++) {
    verts[0] = i;
    std::vector<float> all_geod_distances = dijkstra(m, verts, maxdist);

    // Filter for vertices within maxdist.
    std::vector<float> geod_distances_in_range;
    std::vector<int> geod_indices_in_range;
    for(int j=0; j<all_geod_distances.size(); j++) {
      if(i==j || (all_geod_distances[j] > 0.0 && all_geod_distances[j] < maxdist)) {
        geod_indices_in_range.push_back(j);
        geod_distances_in_range.push_back(all_geod_distances[j]);
      }
    }
    geod_indices.push_back(geod_indices_in_range);
    geod_distances.push_back(geod_distances_in_range);
  }

  std::vector<std::vector<float>> gaussian_weights = gauss_weights(geod_indices, geod_distances, gstd);
  std::vector<float> smoothed_data = spatial_filter(data, geod_indices, gaussian_weights);
  return wrap(smoothed_data);
}

