#if !defined(__NUMERICAL_OPTIMIZATION_H__)
#define __NUMERICAL_OPTIMIZATION_H__
#include <functional>
#include "mt_log.h"
namespace mt {
double NROptimize(double initial, // init branch lengths
                  std::size_t _maxiter,
                  const std::function<double(double, double &, double &)> lnLAnd2Derivs,
                  const double minParam, 
                  const double maxParam);
} // namespace
#endif
