
#include <RcppArmadillo.h>

using namespace Rcpp;


RcppExport SEXP Rpvd_smoothnn(SEXP data_, SEXP num_iter_, SEXP neighborhood_)
{
  try {
    // Convert variables
    std::vector<NumericVector> neighborhood = Rcpp::as<std::vector<NumericVector>>(neighborhood_);
    int num_iter = Rcpp::as<int>(num_iter_);
    NumericVector data(data_);

    // Run smoothing.
    NumericVector smoothed;
    NumericVector source_data;
    float sum, cur_val;
    int neigh_size; // Size of considered neighbors, i.e., only those with non-NA values count.
    for(int iter_idx=0; iter_idx<num_iter; iter_idx++) {
      source_data = i == 0 ? data : smoothed;
      for(int vidx=0; vidx<data.length(); vidx++) {
        if(NumericVector::is_na(source_data[vidx])) {
          continue;
        } else {
          sum = 0.0f;
          neigh_size = 0;
          for(int neigh_idx=0; neigh_idx<neighborhood[vidx].length(); neigh_idx++) {
            cur_val = data[neighborhood[vidx][neigh_idx]];
            if(! NumericVector::is_na(cur_val)) {
              sum += cur_val;
              neigh_size++;
            }
          }
          if(neigh_size > 0) {
            smoothed[vidx] = sum / neigh_size;
          } else {
            smoothed[vidx] = source_data[vidx];
          }

        }
      }
    }
    return wrap(smoothed);

  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

