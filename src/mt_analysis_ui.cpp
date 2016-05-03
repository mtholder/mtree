#include "mt_instance.h"
#include "mt_optimize_branches.h"
#include "mt_likelihood.h"
#include "search.h"
#include "pattern_class.h"
#include "mt_opt_model.h"
#include "mt_ini_options.h"
#include <iostream>

namespace mt {

void doAnalysis(std::ostream * os, MTInstance & instance, const INIBasedSettings & ibs) {
    const enum ProcessActionsEnum action = ibs.action;
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
      double diff = 0.0;
      int howmanyloops = 0;
      do {
        prevlnL = ScoreTree(instance.partMat, instance.tree, instance, false);
        //optimizeModel(instance,.1);
        optimizeModelUsingGolden(instance);
        currlnL = optimizeAllBranchLengths(instance);
        howmanyloops++;
        diff = prevlnL - currlnL;
        std::cerr << "lnL difference for high level opt loop: " << diff << "\n";
      } while(fabs(prevlnL - currlnL) > lnLEpsilon);
      double endlnL = currlnL;
      std::cerr << "lnL before full optimization = " << startlnL << "\n";
      std::cerr << "lnL after full optimization = " << endlnL << "\n";
      *os << endlnL << "\t";
      std::cerr << "Number of high level optimization loops: " << howmanyloops << "\n";
      //instance.tree.write(*os);
    } else if (action == TREE_SEARCH) {
      simpleSPRSearch(instance, 100);
      /*
      //int steps = 10;
      double startL = ScoreTree(instance.partMat, instance.tree, instance, false);
      *os << "Starting likelihood = " << startL << "\n";
      Node * p = instance.tree.GetLeaf(4)->parent->parent;
      Node *par = p->parent;
      Node *subt = removeSubTree(instance, p);
      std::cout << "Got here\n";
      double endL = ScoreTree(instance.partMat, instance.tree, instance, false);
      *os << "Likelihood after subtree removed = " << endL << "\n";
      //assert(p->parent);
      p = insertSubTree(instance, subt, instance.tree.GetRoot()->leftChild->rightSib, par);
      double nextL = ScoreTree(instance.partMat, instance.tree, instance, false);
      *os << "Likelihood after subtree inserted = " << nextL << "\n";
      */
    //  performSearch(instance, steps, instance.tree);
  } else if (action == MISC_TEST) {
    //std::cout << "Made it here\n";
    double result = totalInformativePatternProb(instance);
    std::cerr << "Probability of only informative patterns = " << result << "\n";
    std::cerr << "lnL of only informative patters = " << log(result) << "\n";
  }
}


} //namespace
