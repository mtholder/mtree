#include "mt_tree.h"
#include "mt_optimize_branches.h"
#include "search.h"

namespace mt {

void copyTree(Tree &tree, int &currndnum, const Node * treenode, Node * nodecopy) {
    if(!treenode->IsLeaf()) {
      // if not a leaf function does nothing because leaf node has already between
      // allocated by previous call

      Node *lc, *rs;
      lc = tree.GetNode(currndnum++);
      rs = tree.GetNode(currndnum++);

      copyNode(treenode->leftChild, lc);
      copyNode(treenode->leftChild->rightSib, rs);

      nodecopy->leftChild = lc;
      nodecopy->leftChild->rightSib = rs;

      lc->parent = nodecopy;
      rs->parent = nodecopy;

      // function recursively called on both children
      copyTree(tree, currndnum, treenode->leftChild, lc);
      copyTree(tree, currndnum, treenode->leftChild->rightSib, rs);
    }
}

void Tree::copyTopology(const Tree &other){
  int currNdNum = other.GetNumLeaves();
  Node * r = this->GetNode(currNdNum++);
  this->SetRoot(r);
  copyTree(*this, currNdNum, other.GetRoot(), r);
}

void Node::write(std::ostream & os) const {
    const Node * c = leftChild;
    if (c) {
        os << '(';
        while (c) {
            if (c != leftChild) {
                os << ',';
            }
            c->write(os);
            c = c->rightSib;
        }
        os << ')';
    }
    if (number < UINT_MAX) {
        os << number;
    }
    os << ':' << edgeLen;
}


} //namespace
