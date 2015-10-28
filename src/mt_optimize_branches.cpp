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

// Implementation of Brent's Method for One-Dimensional Parameter Optimization, based on brentGeneric in PLL
template<typename FUN>
val_lnl_t maximizeScoreForBracketed(FUN fun,
                                    val_lnl_t lower,
                                    val_lnl_t mid,
                                    val_lnl_t upper) {
    val_lnl_t best = mid;
    const double lowerBound = 0.0;
    const double upperBound = MAX_BR_LEN;
    const unsigned MAX_NUM_ITER = 100;
    const double TOL = 0.0001; // hard-coded. TEMP should be runtime
    const double BRENT_Z_EPSILON  = 1.e-5;
    const double BRENT_VAR_VALUE = 0.3819660;
    //initialize variables
    bool allConverged = false;
    bool converged = false;
    double e = 0.0;
    double d = 0.0;
    double a = lower.first;
    double b = upper.first;
    double x = mid.first;
    double w = mid.first;
    double v = mid.first;
    double fw = -mid.second; // negate to maximize scpre while the function's logic minimizes
    double fv = -mid.second;
    double fx = -mid.second;
    double u = 0.0;
    for (auto iter = 0U; iter < MAX_NUM_ITER; ++iter) {
        assert(a >= lowerBound && a <= upperBound);
        assert(b >= lowerBound && b <= upperBound);
        assert(x >= lowerBound && x <= upperBound);
        assert(v >= lowerBound && v <= upperBound);
        assert(w >= lowerBound && w <= upperBound);
        double xm = 0.5 * (a + b);
        double tol1 = TOL * fabs(x) + BRENT_Z_EPSILON;
        double tol2 = 2.0 * tol1;
        if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
            best.second = -fx;
            best.first = x;
            converged = true;
            return best;
        } else {
            if (fabs(e) > tol1) {
                const double r = (x - w) * (fx - fv);
                double q = (w - v) * (fx - fw);
                double p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if (q > 0.0) {
                    p = -p;
                }
                q = fabs(q);
                const double etemp = e;
                e = d;
                if((fabs(p) >= fabs(0.5 * q * etemp))
                   || (p <= q * (a-x)) || (p >= q * (b - x))) {
                    d = BRENT_VAR_VALUE * (e = (x >= xm ? a - x : b - x));
                } else {
                    d = p / q;
                    u = x + d;
                    if (u - a < tol2 || b - u < tol2) {
                        d = (xm - x > 0.0 ? fabs(tol1) : -fabs(tol1));
                    }
                }
            } else {
              d = BRENT_VAR_VALUE * (e = (x >= xm ? a - x : b - x));
            }
            u = ((fabs(d) >= tol1) ? (x + d) : (x + (tol1 > 0.0 ? fabs(d) : -fabs(d))));
        }
        assert(converged || (u >= lowerBound && u <= upperBound));
        double fu = -fun(u); // negate to maximize, while function is a minimizer
        if(fu <= fx) {
            if (u >= x) {
                a = x;
            } else {
                b = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if (u < x){
                a = u;
            } else {
                b = u;
            }
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else {
                if(fu <= fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                }
            }
        }
        assert(a >= lowerBound && a <= upperBound);
        assert(b >= lowerBound && b <= upperBound);
        assert(x >= lowerBound && x <= upperBound);
        assert(v >= lowerBound && v <= upperBound);
        assert(w >= lowerBound && w <= upperBound);
        assert(u >= lowerBound && u <= upperBound);
    }
    return best;
}


double maximizeLnLForBrLen(MTInstance &instance, Arc & arc, double prevScore);



double maximizeLnLForBrLen(MTInstance &instance, Arc & arc, double prevScore) {
    std::cerr << "maximizeLnLForBrLen for node terminating with " << arc.toNode->GetNumber() << ' ...\n';
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
            const auto prevLnL = currLnL;
            currLnL = maximizeLnLForBrLen(instance, arc, currLnL);
            const auto thisArcDiff = currLnL - prevLnL;
            assert(fabs(thisArcDiff) >= 0.0);
            std::cerr << "   next Arc optimizeAllBranchLengths tsi = " << tsi << "  currLnL = " << currLnL << " prevLnL = " << prevLnL << '\n';
            arc = poTrav.next();
        } while(arc.toNode);
        const auto thisRoundImprovement = currLnL - beforeThisRound;
        std::cerr << "optimizeAllBranchLengths tsi = " << tsi << "  currLnL = " << currLnL << " thisRoundImprovement = " << thisRoundImprovement << '\n';
        if (thisRoundImprovement < tolForSweep) {
            return currLnL;
        }
    }
    std::cerr << "WARNING: you have hit placeholder code for optimizeAllBranchLengths\n";
    return currLnL;
}


} // namespace mt

