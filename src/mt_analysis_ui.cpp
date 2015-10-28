#include "mt_instance.h"
#include "mt_optimize_branches.h"
#include "mt_likelihood.h"
#include "search.h"
#include "mt_opt_model.h"
#include <iostream>

namespace mt {

void doAnalysis(std::ostream * os, MTInstance & instance, enum ProcessActionsEnum action) {
    if (action == SCORE_ACTION) {
        const double lnL = ScoreTree(instance.partMat, instance.tree, instance);
        if (os) {
            *os << "lnL = " << lnL << "\n";
        }
    } else if (action == OPTIMIZE_BR_LEN) {
        const double beforelnL = ScoreTree(instance.partMat, instance.tree, instance);
        if (os) {
            *os << "lnL before brlen optimization = " << beforelnL << "\n";
        }
        const double afterLnL = optimizeAllBranchLengths(instance);
        const double diff = afterLnL - beforelnL;
        if (os) {
            *os << "lnL before brlen optimization = " << beforelnL << "\n";
            *os << "lnL after brlen optimization = " << afterLnL << "\n";
            *os << "Change in lnL = " << diff << "\n";
            instance.tree.write(*os);
        }

    } else if (action == OPTIMIZE_PARS) {
      double startlnL = ScoreTree(instance.partMat, instance.tree, instance);
      optimizeModel(instance,.1);
      double endlnL = ScoreTree(instance.partMat, instance.tree, instance);
      *os << "lnL before param optimization: " << startlnL << "\n";
      *os << "lnL after param optimization: " << endlnL << "\n";
    } else if (action == FULL_OPTIMIZE) {
      double prevlnL = ScoreTree(instance.partMat, instance.tree, instance);
      double startlnL = prevlnL;
      double currlnL = ScoreTree(instance.partMat, instance.tree, instance);
      double lnLEpsilon = 0.1;
      do {
        prevlnL = ScoreTree(instance.partMat, instance.tree, instance);
        optimizeModel(instance,.1);
        currlnL = optimizeAllBranchLengths(instance);
      } while(fabs(prevlnL - currlnL) > lnLEpsilon);
      double endlnL = currlnL;
      *os << "lnL before full optimization = " << startlnL << "\n";
      *os << "lnL after full optimization = " << endlnL << "\n";
    } else if (action == TREE_SEARCH) {
      //int steps = 10;
      double startL = ScoreTree(instance.partMat, instance.tree, instance);
      *os << "Starting likelihood = " << startL << "\n";
      Node * p = instance.tree.GetLeaf(4)->parent->parent;
      mtreeTestSPR (instance, p, 2, startL);
      double endL = ScoreTree(instance.partMat, instance.tree, instance);
      *os << "Likelihood after subtree removed = " << endL << "\n";
    //  performSearch(instance, steps, instance.tree);
    }
}


} //namespace
