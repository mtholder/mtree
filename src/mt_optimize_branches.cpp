#include "mt_instance.h"
#include "mt_optimize_branches.h"
#include "mt_likelihood.h"
#include "mt_tree_traversal.h"
namespace mt {
double maximizeScoreForBrLen(MTInstance &instance, Arc & arc, double prevScore);


double maximizeScoreForBrLen(MTInstance &instance, Arc & arc, double prevScore) {
    auto brLenScorer = [&] (double nu) {
        const double prev = arc.GetEdgeLen();
        arc.SetEdgeLen(nu);
        double lnL = ScoreTree(instance.partMat, instance.tree, instance);
        arc.SetEdgeLen(prev);
        return lnL;
    };

    const double tiny = 0.001;
    const double tinyBrLenLnL = brLenScorer(tiny);
    const double mid = 0.05;
    const double midBrLenLnL = brLenScorer(mid);
    const double large = 1.0;
    const double largeBrLenLnL = brLenScorer(large);
    double optBrLen = mid;
    double optBrLenLnL = midBrLenLnL;

    std::cerr << "   nu = " << tiny << " ===> lnL = " << tinyBrLenLnL << '\n';
    std::cerr << "   nu = " << mid << " ===> lnL = " << midBrLenLnL << '\n';
    std::cerr << "   nu = " << large << " ===> lnL = " << largeBrLenLnL << '\n';
    arc.SetEdgeLen(optBrLen);
    return optBrLenLnL;
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
            currLnL = maximizeScoreForBrLen(instance, arc, currLnL);
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

