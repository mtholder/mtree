//Tree Search Algorithm from PLL - searchAlgo.c + topologies.c

#if !defined(__SEARCH_H__)
#define __SEARCH_H__
#include "mt_tree.h"
#include "mt_instance.h"
#include "mt_data.h"

namespace mt {

class searchInfo {
  double bestLnL;
  Tree bestTree;
  bool converged;
  searchInfo(MTInstance &instance)
    :bestLnL(MT_UNLIKELY),
    converged(false),
    bestTree(instance.tree.GetNumNodes(), instance.tree.GetNumLeaves()) {
      initBestTree(instance);
    }
  void initBestTree(MTInstance &instance);
};

/* funcs */

//void searchInfo::initBestTree(MTInstance &instance);
void copyNode(Node *r, Node *s);
void mtreeTestSPR(MTInstance &instance, Node * p, int maxtrav, double bestLnL);
Node * removeSubTree(MTInstance &instance, Node * p);
Node * insertSubTree(MTInstance &instance, Node *p, Node *q, Node *s);

} //namespace mt

#endif
