#include "mt_instance.h"
#include "mt_optimize_branches.h"
#include "mt_likelihood.h"
#include "mt_tree_traversal.h"
namespace mt {

typedef std::pair<double, double> val_lnl_t;
template<typename FUN>
val_lnl_t maximizeScoreForSmallBrLen(FUN fun, val_lnl_t);
template<typename FUN>
val_lnl_t maximizeScoreForLargeBrLen(FUN fun, val_lnl_t);
template<typename FUN>
val_lnl_t maximizeScoreForBracketed(FUN fun, val_lnl_t lower, val_lnl_t mid, val_lnl_t upper);


template<typename FUN>
val_lnl_t maximizeScoreForSmallBrLen(FUN fun, val_lnl_t prev) {
    val_lnl_t best = prev;
    return best;
}

template<typename FUN>
val_lnl_t maximizeScoreForLargeBrLen(FUN fun, val_lnl_t prev) {
    val_lnl_t best = prev;
    return best;
}

template<typename FUN>
val_lnl_t maximizeScoreForBracketed(FUN fun, val_lnl_t lower, val_lnl_t mid, val_lnl_t upper) {
    val_lnl_t best = mid;
    return best;
}


double maximizeLnLForBrLen(MTInstance &instance, Arc & arc, double prevScore);



double maximizeLnLForBrLen(MTInstance &instance, Arc & arc, double prevScore) {
    auto brLenScorer = [&] (double nu) {
        const double prev = arc.GetEdgeLen();
        arc.SetEdgeLen(nu);
        double lnL = ScoreTree(instance.partMat, instance.tree, instance);
        arc.SetEdgeLen(prev);
        return lnL;
    };
    val_lnl_t soln{arc.GetEdgeLen(), prevScore};
    const val_lnl_t tiny{0.001, brLenScorer(0.001)};
    const val_lnl_t small{0.01, brLenScorer(0.01)};
    std::cerr << "   nu = " << tiny.first << " ===> lnL = " << tiny.second << '\n';
    std::cerr << "   nu = " << small.first << " ===> lnL = " << small.second << '\n';
    if (tiny.second > small.second) {
        soln = maximizeScoreForSmallBrLen(brLenScorer, tiny);
    } else {
        const val_lnl_t mid{0.05, brLenScorer(0.05)};
        std::cerr << "   nu = " << mid.first << " ===> lnL = " << mid.second << '\n';
        if (small.second > mid.second) {
            soln = maximizeScoreForBracketed(brLenScorer, tiny, small, mid);
        } else {
            const val_lnl_t large{1.00, brLenScorer(1.00)};
            std::cerr << "   nu = " << large.first << " ===> lnL = " << large.second << '\n';
            if (mid.second >= large.second) {
                soln = maximizeScoreForBracketed(brLenScorer, small, mid, large);
            } else {
                const val_lnl_t huge{100.00, brLenScorer(100.00)};
                std::cerr << "   nu = " << huge.first << " ===> lnL = " << huge.second << '\n';
                if (large.second > huge.second) {
                    soln = maximizeScoreForBracketed(brLenScorer, mid, large, huge);
                } else {
                    soln = maximizeScoreForLargeBrLen(brLenScorer, huge);
                }
            }
        }
    }
    arc.SetEdgeLen(soln.first);
    return soln.second;
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
            currLnL = maximizeLnLForBrLen(instance, arc, currLnL);
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

