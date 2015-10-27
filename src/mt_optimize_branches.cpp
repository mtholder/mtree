#include "mt_instance.h"
#include "mt_optimize_branches.h"
#include "mt_likelihood.h"
#include "mt_tree_traversal.h"
namespace mt {

typedef std::pair<double, double> val_lnl_t;
template<typename FUN>
val_lnl_t maximizeScoreForSmallBrLen(FUN fun, val_lnl_t better, val_lnl_t upper);
template<typename FUN>
val_lnl_t maximizeScoreForLargeBrLen(FUN fun, val_lnl_t better, val_lnl_t lower);
template<typename FUN>
val_lnl_t maximizeScoreForBracketed(FUN fun, val_lnl_t lower, val_lnl_t mid, val_lnl_t upper);

const double SAME_LN_L_TOL = 1e-04;
template<typename FUN>
val_lnl_t maximizeScoreForSmallBrLen(FUN fun, val_lnl_t curr, val_lnl_t upper) {
    val_lnl_t best = curr;
    val_lnl_t zeroBound{0.0, fun(0.0)};
    if (zeroBound.second > curr.second) {
        double nextBrLen = curr.first/2.0;
        val_lnl_t next{nextBrLen, fun(nextBrLen)};
        if (next.second > zeroBound.second) {
            return maximizeScoreForBracketed(fun, zeroBound, next, curr);
        }
        while (zeroBound.second - next.second > SAME_LN_L_TOL) {
            curr = next;
            next.first = next.first / 2.0;
            next.second = fun(next.first);
            if (next.second > zeroBound.second) {
                return maximizeScoreForBracketed(fun, zeroBound, next, curr);
            }
        }
        best = zeroBound;
    } else {
        return maximizeScoreForBracketed(fun, zeroBound, curr, upper);
    }
    return best;
}

const double MAX_BR_LEN = 300.0;
template<typename FUN>
val_lnl_t maximizeScoreForLargeBrLen(FUN fun, val_lnl_t curr, val_lnl_t lower) {
    val_lnl_t best = curr;
    val_lnl_t upperBound{MAX_BR_LEN, fun(MAX_BR_LEN)};
    if (upperBound.second > curr.second) {
        double nextBrLen = (curr.first + upperBound.first)/2.0;
        val_lnl_t next{nextBrLen, fun(nextBrLen)};
        if (next.second > upperBound.second) {
            return maximizeScoreForBracketed(fun, lower, next, upperBound);
        }
        while (upperBound.second - next.second > SAME_LN_L_TOL) {
            curr = next;
            next.first = (next.first + upperBound.first)/2.0;
            next.second = fun(next.first);
            if (next.second > upperBound.second) {
                return maximizeScoreForBracketed(fun, curr, next, upperBound);
            }
        }
        best = upperBound;
    } else {
        return maximizeScoreForBracketed(fun, lower, curr, upperBound);
    }
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
        soln = maximizeScoreForSmallBrLen(brLenScorer, tiny, small);
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
                    soln = maximizeScoreForLargeBrLen(brLenScorer, huge, large);
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

