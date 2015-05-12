#include "mt_numerical_optimization.h"
#include "mt_util.h"
#include <cmath>
using namespace std;
namespace mt {
/*  topLevelMakenewz modified to make a generic, one-parameter optimizer....
   the function below actually implements the iterative Newton-Raphson procedure.
   It is particularly messy and hard to read because for the case of per-partition branch length 
   estimates it needs to keep track of whetehr the Newton Raphson procedure has 
   converged for each partition individually. 
   The rationale for doing it like this is also provided in:
   A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009,
*/
double NROptimize(double z0, // init branch lengths
                  std::size_t _maxiter,
                  const std::function<double(double, double &, double &)> lnLAnd2Derivs,
                  const double minParam, 
                  const double maxParam) {
      // initialize loop convergence variables etc. 
      // maxiter is the maximum number of NR iterations we are going to do before giving up 
    double z = z0;
    std::size_t maxiter = _maxiter;
    bool curvatOK= false;
    double zstep;
    double zprev = z;
    for (;;) {
        if (curvatOK) {
            curvatOK = false;
            zprev = z;
            zstep = (1.0 - maxParam) * z + minParam;
        }
        enforce_bound(z, minParam, maxParam);
        double dlnLdlz, d2lnLdlz2;
        lnLAnd2Derivs(z, dlnLdlz, d2lnLdlz2);
            // do a NR step, if we are on the correct side of the maximum that's okay, otherwise 
            // shorten branch 
        if(!curvatOK) {
            if ((d2lnLdlz2 >= 0.0) && (z < maxParam)) {
                // Bad curvature, shorten branch
                z = 0.37*z + 0.63;
                zprev = z;
            } else {
                curvatOK = true;
            }
        }
            // do the standard NR step to obrain the next value, depending on the state for eahc partition
        if (curvatOK) {
            const double pb = 0.25 * zprev + 0.75;
            if (d2lnLdlz2 < 0.0) {
                double tantmp = -dlnLdlz / d2lnLdlz2;
                if (tantmp < 100) {
                    z *= exp(tantmp);
                }
                enforce_min_bound(z, minParam);
                enforce_max_bound(z, pb);
            } else {
                z = pb;
            }
        }
        enforce_max_bound(z, maxParam);
        if (std::abs(z - zprev) <= zstep) {
            break;
        }
        --maxiter;
           // We should make a more informed decision here, based on the log like improvement
        if (maxiter == 0) {
            z = z0;
            break;
        }
    }
    return z;
}

} //namespace