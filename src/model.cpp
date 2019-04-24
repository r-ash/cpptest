#include "model.h"

// Notes before putting down there are some issues here with compilation. The problems are arising because entrantPrev may be null in the R code
// overloading methods here doesn't seem to work but also trying to use boost::optional<> types doesn't work as it isn't an SEXP type so Rcpp
// cannot convert it for being alled from R. Using Rcpp::Nullable type also doesn't work because this can't wrap an std::vector.