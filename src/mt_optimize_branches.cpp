#include "mt_optimize_branches.h"
#include "mt_numerical_optimization.h"
#include "mt_instance.h"
#include <cassert>
namespace mt {

class EdgeLenOptAdaptor {
    public:
    EdgeLenOptAdaptor(MTInstance &, Arc & ) {

    }
    double operator()(double z, double & dlnLdlz, double & d2lnLdlz2) {
        dlnLdlz = 0.0;
        d2lnLdlz2 = -1.0;
        return 0.0;
    }
};
double optimizeAllLengthsForOneEdge(MTInstance & mtInstance, Arc & edge) {
    const auto z0 = edge.GetEdgeLen();
    const std::size_t maxiter = 21;
    EdgeLenOptAdaptor adaptor(mtInstance, edge);
    auto z = NROptimize(z0, maxiter, adaptor, PLL_ZMIN, PLL_ZMAX);
    return z;
}

/*
    auto coreLogZ = std::log(z);
    SetExecuteMaskByBranchLengthPart(instance, curvatOK, true);
    StoreExecuteMaskInTraversalDescriptor(instance);
    StoreValuesInTraversalDescriptor(instance, coreLogZ);
    double dlnLdlz, d2lnLdlz2;
        // sequential part, if this is the first newton-raphson implementation,
        // do the precomputations as well, otherwise just execute the computation
        // of the derivatives
    if (firstIteration) {
        makenewzIterative(instance);
        firstIteration = false;
    }
    execCore(instance, dlnLdlz, d2lnLdlz2);
*/
} // namespace