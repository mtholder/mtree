//Tree Search Algorithm from PLL - searchAlgo.c + topologies.c

#if !defined(__SEARCH_H__)
#define __SEARCH_H__
#include "mt_tree.h"

class BestTree {
  double LnL;
  Tree tree;
  BestTree(double likelihood, Tree tr)
    :LnL(likelihood),
    Tree(tr) {
    }
}

class BTList {
  BestTree *list;
  int n;
  BTList(BestTree *list, int n)
    :BestTree(list),
    int(n) {
    }
}

class rearrangeStep {
    public:
      enum rearrangeType {
          SPR_MODE = 0,
          NNI_MODE = 1
          // TBR_MODE = 2
      };
      double  likelihood;
      union   rearrangeType {
        class SPR {
          Node * removeNode;
          Node * insertNode;
          double zqr[NUM_BRANCHES];
          };
        class NNI {
            Node * originNode;
            int    swapType;
          };
      /*class TBR {
      }*/
      };
};

class nni_Move {
        MTInstance * tr;
        Node * p;
        int nniType;
        double likelihood;
        double deltaLH;
};

//Funcs

/* rearrange functions (NNI and SPR) */
void performSearch(MTInstance &instance, int steps, Node *p);
void searchStep(MTInstance &instance);
int bestSPR(MTInstance &instance, Node *p, int mintrav, int maxtrav);
int performNNI(MTInstance &instance);
