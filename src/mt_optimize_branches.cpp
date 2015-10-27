#include "mt_instance.h"
#include "mt_optimize_branches.h"
#include "mt_likelihood.h"
#include "mt_tree_traversal.h"
namespace mt {
double optimizeSingleBranchLength(MTInstance & instance, Arc & arc, double prevLnL);

double optimizeSingleBranchLength(MTInstance & instance, Arc & arc, double ) {
    auto & tree = instance.tree;
    arc.SetEdgeLen(0.0);
    double zeroBrLenLnL = ScoreTree(instance.partMat, tree, instance);
    arc.SetEdgeLen(100);
    double saturatedBrLenLnL = ScoreTree(instance.partMat, tree, instance);
    arc.SetEdgeLen(0.05);
    double ptZeroFiveBrLenLnL = ScoreTree(instance.partMat, tree, instance);
    std::cerr << "   nu = 0.000 ===> lnL = " << zeroBrLenLnL << '\n';
    std::cerr << "   nu = 0.050 ===> lnL = " << ptZeroFiveBrLenLnL << '\n';
    std::cerr << "   nu = 100.0 ===> lnL = " << saturatedBrLenLnL << '\n';
    return ptZeroFiveBrLenLnL;
}

double optimizeAllBranchLengths(MTInstance &instance) {
    const unsigned maxNumTreeSweeps = 10; // TEMP! should be runtime
    const double tolForSweep = 0.0001; // TEMP! should be runtime
    auto & tree = instance.tree;
    auto rootPtr = tree.GetRoot();
    assert(rootPtr);
    const double beforeOpt = ScoreTree(instance.partMat, tree, instance);
    double currLnL = beforeOpt;
    for (auto tsi = 0U; tsi < maxNumTreeSweeps; ++tsi) {
        const auto beforeThisRound = currLnL;
        PostorderForNodeIterator poTrav = postorder(rootPtr);
        Arc arc = poTrav.get();
        do {
            currLnL = optimizeSingleBranchLength(instance, arc, currLnL);
            arc = poTrav.next();
        } while(arc.toNode);
        const auto thisRoundImprovement = currLnL - beforeThisRound;
        if (thisRoundImprovement < tolForSweep) {
            return currLnL;
        }
    }
    std::cerr << "WARNING: you have hit placeholder code for optimizeAllBranchLengths\n";
    return currLnL;
}


} // namespace mt

