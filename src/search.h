//Tree Search Algorithm from PLL - searchAlgo.c + topologies.c

#if !defined(__SEARCH_H__)
#define __SEARCH_H__
#include "mt_tree.h"
#include "mt_instance.h"

namespace mt {

/*class BestTree {
  double LnL;
  Tree tree;
  BestTree(double likelihood, Tree tr)
    :LnL(likelihood),
    Tree(tr) {
    }
};
*/
/* rearrange functions (NNI and SPR) */
void performSearch(MTInstance &instance, int steps, Node *p);
void searchStep(MTInstance &instance);
int bestSPR(MTInstance &instance, Node *p, int mintrav, int maxtrav);
int performNNI(MTInstance &instance);
void mtreeTestSPR(MTInstance &instance, Node * p, int maxtrav, double bestLnL);
Node * removeSubtree(MTInstance &instance, Node * p);

} //namespace mt

#endif
