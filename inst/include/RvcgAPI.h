#include <Rcpp.h>

using namespace Rcpp;


#ifndef _RVCG_API_H
#define _RVCG_API_H

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

//#include <xts.h>		// also includes R.h, Rinternals.h, Rdefines.h

#include <Rconfig.h>
#include <R_ext/Rdynload.h>

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

#ifdef __cplusplus
  extern "C" {
#endif


SEXP attribute_hidden Rdijkstra(SEXP vb_, SEXP it_, SEXP verts_, SEXP maxdist_) {
  static SEXP(*fun)(SEXP,SEXP,SEXP,SEXP) = NULL;
  if (fun == NULL)
    fun = (SEXP(*)(SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("Rvcg","Rdijkstra");
  return fun(x, k, pad);
}

#ifdef __cplusplus
}
#endif

#endif /* !_RVCG_API_H */

