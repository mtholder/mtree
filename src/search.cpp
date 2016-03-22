#include "search.h"
#include "mt_tree.h"
#include "mt_likelihood.h"

#include <assert.h>
#include <vector>

// Implementation of Tree Topology Search Algorithm ported from PLL (searchAlgo.c and topologies.c)
using namespace std;
namespace mt {

// copies basic info of node l into node k
void copyNode(Node * l, Node * k) {
    //std::cout << "Copying node " << l->GetNumber() << "\n";
    k->SetNumber(l->GetNumber());
    k->SetEdgeLen(l->GetEdgeLen());
}

/*
// frees tree data structure
// recursively deletes nodes in tree

void freeTree(Node * nd){

  if(nd->IsLeaf()) {
    delete nd;
    //nd = NULL;

  } else {

    Node * lc = nd->leftChild;
    Node * rs = lc->rightSib;
    delete nd;
    //nd = NULL;

    // recursive call
    freeTree(lc);
    freeTree(rs);
  }
}
*/

// recursively copies tree structure and basic node data
// from nodes starting from treenode to nodes starting at nodecopy (roots)

// is recursion best for this?



// allocate memory for best tree topology to be saved
// and copy from current tree
void searchInfo::initBestTree(MTInstance &instance){

  Tree * scratch = new Tree(instance.tree);
}
  /*Node * root = instance.tree.GetRoot();
  copyNode(root, r);

  hardCopyTree(root, r);
  */


// remove and return subtree located at node p, collapse branch on other tree
// should optimize new branch length between p->parent->parent and p's sibling after removal
// N.B. need to pass parent of p to a placeholder before calling this function
// otherwise the removed node will be lost!
Node * removeSubTree(MTInstance & instance, Node * p) {
  //cout << "removeSubTree called\n";
  assert (p);
  assert (p->parent);
  assert (p->parent->parent); // parent must have a parent b/c root is constrained
  //cout << p->parent->parent << "\n";
  // assert(IsInTree(instance, p));
  //if (!p.IsLeaf()) {} should not matter if it is leaf or internal node
  if (p->parent->IsLeftChild()) {
    if (p->IsLeftChild()) {
      //cout <<"left child of a left child\n";
      p->parent->parent->leftChild = p->rightSib;
      p->rightSib->parent = p->parent->parent;
      p->rightSib->rightSib = p->parent->rightSib;
      p->rightSib = nullptr;
    } else {
      //cout << "right child of a left child\n";
      p->parent->parent->leftChild = p->parent->leftChild;
      p->parent->leftChild->parent = p->parent->parent;
      p->parent->leftChild->rightSib = p->parent->rightSib;
    }
  } else {
    if (p->IsLeftChild()) {
      //cout << "left child of a right child\n";
      p->parent->parent->leftChild->rightSib = p->rightSib;
      p->rightSib->parent = p->parent->parent;
      p->rightSib = nullptr;
    } else {
      //cout << "right child of a right child\n";
      p->parent->parent->leftChild->rightSib = p->parent->leftChild;
      p->parent->leftChild->parent = p->parent->parent;
      p->parent->leftChild->rightSib = nullptr;
    }
  }

  p->parent = nullptr;
  return p;
}

// inserts subtree p as rightSib to node q on main tree
// inserts saved node s as parent to q and p
// returns root to tree
Node * insertSubTree(MTInstance &instance, Node * p, Node * q, Node *s) {
  assert(q->parent); // for now only dealing with rooted trees
  assert(s);
  bool lc = q->IsLeftChild();
  //std::cout << lc << "\n";
  // these pointers don't depend on whether q is a right or left child
  p->parent = s;
  s->leftChild = q;
  s->parent = q->parent;

  if (lc){
    s->rightSib = q->rightSib;
    q->parent->leftChild = s;
  } else {
    q->parent->leftChild->rightSib = s;
    s->rightSib = nullptr;
  }

  q->rightSib = p;
  q->parent = s;

  assert(q->parent);
  assert(s->parent);
  assert(q->rightSib);


  return instance.tree.GetRoot();
}

void simpleSPRSearch(MTInstance &instance, int maxloops) {
  searchInfo sInfo(instance);
  instance.tree.TreeDebug();
  sInfo.bestLnL = ScoreTree(instance.partMat, instance.tree, instance, false);
  std::cout << "Starting ln likelihood = " << sInfo.bestLnL << "\n";
  bool changed = false;
  int step = 0;
  int location = 0;
  while(step++ < maxloops) {

    std::cout << "Now on loop # " << step << "\n";

    for(int i = 0; i < instance.tree.GetNumNodes(); i++) {
      std::cout << "Trying node " << i << "\n";
      if (!instance.tree.GetNode(i)->parent) continue;
      if (!instance.tree.GetNode(i)->parent->parent) continue;
      Node * snipNode = instance.tree.GetNode(i);

      if (snipNode->IsLeftChild()) {
        location = snipNode->rightSib->GetNumber();
      } else {
        location = snipNode->parent->leftChild->GetNumber();
      }

      Node * temp = snipNode->parent;
      Node * subt = removeSubTree(instance, snipNode);

      double templnl = ScoreTree(instance.partMat, instance.tree, instance, false);
      std::cout << "Lnl of pruned tree: " << templnl << "\n";
      //instance.tree.TreeDebug();

      changed = false;

      // try insertions at all nodes not in subtree
      for(int j = 0; j < instance.tree.GetNumNodes(); j++) {
        if (instance.tree.GetNode(j)->parent == nullptr ||
            instance.tree.GetNode(j)->parent == instance.tree.GetRoot()) {
              //std::cout << "Not valid insertion point\n";
              continue;
            } else {
              if(instance.tree.isNodeConnected(instance.tree.GetRoot(), j)) {
                std::cout << "Trying insertion at node " << j << "\n";
                Node * newroot = insertSubTree(instance, subt, instance.tree.GetNode(j), temp);
                instance.tree.TreeDebug();
                // only scores new lnl if subtree has been reinserted
                //std::cout << "Scoring new tree\n";
                double newlnl = ScoreTree(instance.partMat, instance.tree, instance, false);
                if (newlnl > sInfo.bestLnL) {
                    sInfo.bestLnL = newlnl;
                    std::cout << "New ln likelihood = " << sInfo.bestLnL << "\n";
                    //sInfo.bestTree.copyTopology(instance.tree);
                    changed = true;
                    break;
                  } else {
                    //temp = snipNode->parent; // is this necessary?
                    subt = removeSubTree(instance, snipNode);
                  }
              }
          }
        }
      // if likelihood has not improved, reinsert to original location
      if (!changed) Node * resetroot = insertSubTree(instance, subt, instance.tree.GetNode(location), temp);
      std::cout << "Got to end of insertion for loop\n";
      //std::cout << "Debugging bestTree\n";
      //sInfo.bestTree.TreeDebug();
      //std::cout << "Debugging instance tree\n";
      //instance.tree.TreeDebug();
      //instance.tree.copyTopology(sInfo.bestTree);
      instance.tree.TreeDebug();
      //sInfo.bestTree.TreeDebug();
    }
    // break out of while loop if there is no improvement in lnl
    if(!changed) break;
  }
  std::cout << "End ln likelihood = " << sInfo.bestLnL << "\n";
  std::cout << "Total steps: " << step << "\n";
}

/*
// Tries SPR moves for a subtree rooted at node p up to maxtrav nodes away.
// Derived from pllTestSPR.
void mtreeTestSPR (MTInstance &instance,
                   Node * p,
                   int maxtrav,
                   double bestLnL
                   ) {
  assert(false);
  //if (1) {
    //double bl = p->edgeLen;
    //Node * n = p->parent;
    //Node * fromnode = p->parent->parent;
    //Node * subtree =
                     removeSubTree(instance, p);
  }
}
*/

} //namespace mt
