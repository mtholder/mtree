#include "mt_instance.h"
#include "mt_optimize_branches.h"
#include "mt_likelihood.h"
namespace mt {

double optimizeAllBranchLengths(MTInstance &instance) {
    std::cerr << "WARNING: you have hit placeholder code for optimizeAllBranchLengths\n";
    return ScoreTree(instance.partMat, instance.tree, instance);
}


} // namespace mt

