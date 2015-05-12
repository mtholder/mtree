#include "search.h"
#include "mt_tree.h"

#include <assert.h>
#include <vector>

// Implementation of Tree Topology Search Algorithm ported from PLL (searchAlgo.c and topologies.c)

//Currently does not use any branch length optimizations
using namespace std;
namespace mt {


Node * removeSubTree(MTInstance & ,//instance,
                     Node * p) {
  cout << "removeSubTree called\n";
  cout << p->parent->parent << "\n";
    if (1) /*(!p.IsLeaf())*/ {
      // assert(IsInTree(instance, p));
      if (p->parent == p->parent->parent->leftChild) {
        if (p == p->parent->leftChild) {
          cout <<"left child of a left child\n";
          p->parent->parent->leftChild = p->parent->rightSib;
          p->parent->rightSib->parent = p->parent->parent;
        } else {
          cout << "right child of a left child\n";
          p->parent->parent->leftChild = p->parent->leftChild;
          p->parent->leftChild->parent = p->parent->parent;
        }
      } else {
        if (p == p->parent->leftChild) {
          cout << "left child of a right child\n";
          p->parent->parent->rightSib = p->parent->rightSib;
          p->parent->rightSib->parent = p->parent->parent;
        } else {
          cout << "right child of a right child\n";
          p->parent->parent->rightSib = p->parent->leftChild;
          p->parent->leftChild->parent = p->parent->parent;
        }
      }
      p->parent = NULL;
    }
    return p;
}

    // Tries SPR moves for a subtree rooted at node p up to maxtrav nodes away.
    // Derived from pllTestSPR.
void mtreeTestSPR (MTInstance &instance, Node * p, int maxtrav, double bestLnL) {
  if (1)/*(!IsLeaf(p))*/ {
    double bl = p->edgeLen;
    Node * n = p->parent;
    Node * fromnode = p->parent->parent;
    Node * subtree = removeSubTree(instance, p);
    return;
  }
}

} //namespace mt
