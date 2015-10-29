#include "mt_instance.h"
#include "mt_optimize_branches.h"
#include "mt_likelihood.h"
#include "search.h"
#include "mt_opt_model.h"
#include <iostream>

namespace mt {

void doAnalysis(std::ostream * os, MTInstance & instance, enum ProcessActionsEnum action) {
    if (action == SCORE_ACTION) {
        const double lnL = ScoreTree(instance.partMat, instance.tree, instance, true);
        if (os) {
            *os << "lnL = " << lnL << "\n";
        }
    } else if (action == OPTIMIZE_BR_LEN) {
        const double beforelnL = ScoreTree(instance.partMat, instance.tree, instance, true);
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
      double startlnL = ScoreTree(instance.partMat, instance.tree, instance, true);
      //optimizeModel(instance,.1);
      optimizeModelUsingGolden(instance);
      double endlnL = ScoreTree(instance.partMat, instance.tree, instance, false);
      *os << "lnL before param optimization: " << startlnL << "\n";
      *os << "lnL after param optimization: " << endlnL << "\n";
    } else if (action == FULL_OPTIMIZE) {
      double prevlnL = ScoreTree(instance.partMat, instance.tree, instance, true);
      double startlnL = prevlnL;
      double currlnL = ScoreTree(instance.partMat, instance.tree, instance, false);
      double lnLEpsilon = 0.1;
      int howmanyloops = 0;
      do {
        prevlnL = ScoreTree(instance.partMat, instance.tree, instance, false);
        //optimizeModel(instance,.1);
        optimizeModelUsingGolden(instance);
        currlnL = optimizeAllBranchLengths(instance);
        howmanyloops++;
      } while(fabs(prevlnL - currlnL) > lnLEpsilon);
      double endlnL = currlnL;
      *os << "lnL before full optimization = " << startlnL << "\n";
      *os << "lnL after full optimization = " << endlnL << "\n";
      *os << "Number of high level optimization loops: " << howmanyloops << "\n";
      instance.tree.write(*os);
    } else if (action == TREE_SEARCH) {
      //int steps = 10;
      double startL = ScoreTree(instance.partMat, instance.tree, instance, false);
      *os << "Starting likelihood = " << startL << "\n";
      Node * p = instance.tree.GetLeaf(4)->parent->parent;
      mtreeTestSPR (instance, p, 2, startL);
      double endL = ScoreTree(instance.partMat, instance.tree, instance, false);
      *os << "Likelihood after subtree removed = " << endL << "\n";
    //  performSearch(instance, steps, instance.tree);
    }
}


} //namespace
