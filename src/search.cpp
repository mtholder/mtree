#include "search.h"
#include "mt_tree.h"

#include <assert.h>
#include <vector>

// Implementation of Tree Topology Search Algorithm ported from PLL (searchAlgo.c and topologies.c)
using namespace std;
namespace mt {

// copies basic info of node l into node k
void copyNode(Node * l, Node * k) {
    k->SetNumber(l->GetNumber());
    k->SetEdgeLen(l->GetEdgeLen());
}


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

// recursively copies tree structure and basic node data
// from nodes starting from treenode to nodes starting at nodecopy (roots)

// is recursion best for this?

void hardCopyTree(Node * treenode, Node * nodecopy) {
  if(!treenode->IsLeaf()) {
    // if not a leaf function does nothing because leaf node has already between
    // allocated by previous call

    Node *lc, *rs;
    lc = new Node;
    rs = new Node;

    copyNode(treenode->leftChild, lc);
    copyNode(treenode->leftChild->rightSib, rs);

    nodecopy->leftChild = lc;
    nodecopy->leftChild->rightSib = rs;

    lc->parent = nodecopy;
    rs->parent = nodecopy;

    // function recursively called on both children
    hardCopyTree(treenode->leftChild, lc);
    hardCopyTree(treenode->leftChild->rightSib, rs);
  }
}

// allocate memory for best tree topology to be saved
// and copy from current tree
void searchInfo::initBestTree(MTInstance &instance){

  Node * r = new Node;
  Node * root = instance.tree.GetRoot();
  copyNode(root, r);

  hardCopyTree(root, r);
}

// remove and return subtree located at node p, collapse branch on other tree
// should optimize new branch length between p->parent->parent and p's sibling after removal
// N.B. need to pass parent of p to a placeholder before calling this function
// otherwise the removed node will be lost!
Node * removeSubTree(MTInstance & instance, Node * p) {
  //cout << "removeSubTree called\n";
  assert (p);
  assert (p->parent);
  assert (p->parent->parent); // parent must have a parent b/c root is constrained
  cout << p->parent->parent << "\n";
  // assert(IsInTree(instance, p));
  //if (!p.IsLeaf()) {} should not matter if it is leaf or internal node
  if (p->parent == p->parent->parent->leftChild) {
    if (p == p->parent->leftChild) {
      //cout <<"left child of a left child\n";
      p->parent->parent->leftChild = p->rightSib;
      p->rightSib->parent = p->parent->parent;
      p->rightSib = nullptr;
    } else {
      //cout << "right child of a left child\n";
      p->parent->parent->leftChild = p->parent->leftChild;
      p->parent->leftChild->parent = p->parent->parent;
      p->parent->leftChild->rightSib = p->parent->parent->leftChild->rightSib;
    }
  } else {
    if (p == p->parent->leftChild) {
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

// inserts subtree p below node q on main tree
// sets p as rightSib of q (arbitrary)
// uses s placeholder for new node parent to p
// returns root to tree
Node * insertSubTree(MTInstance &instance, Node * p, Node * q, Node *s) {
  assert(q->parent); // for now only dealing with rooted trees

  // these pointers don't depend on whether q is a right or left child
  p->parent = s;
  s->leftChild = q;
  s->parent = q->parent;

  if (q->IsLeftChild()){
    s->rightSib = q->rightSib;
    q->parent->leftChild = s;
  } else {
    q->parent->leftChild->rightSib = s;
  }

  q->rightSib = p;
  q->parent = s;
   //
}

// Tries SPR moves for a subtree rooted at node p up to maxtrav nodes away.
// Derived from pllTestSPR.
void mtreeTestSPR (MTInstance &instance,
                   Node * p,
                   int maxtrav,
                   double bestLnL
                   ) {
  assert(false);
  if (1)/*(!IsLeaf(p))*/ {
    //double bl = p->edgeLen;
    //Node * n = p->parent;
    //Node * fromnode = p->parent->parent;
    //Node * subtree =
                     removeSubTree(instance, p);
  }
}

} //namespace mt
