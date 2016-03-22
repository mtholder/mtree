#include "mt_tree.h"
#include "mt_tree_traversal.h"
#include "mt_optimize_branches.h"
#include "search.h"

namespace mt {

void copyTree(Tree &tree, const Node * treenode, Node * nodecopy) {
    if(!treenode->IsLeaf()) {
      // if a leaf function does nothing because leaf node has already between
      // allocated by previous call

      //std::cout << "Beginning of copyTree\n";
      Node *lc, *rs;
      lc = tree.GetNode(treenode->leftChild->GetNumber());
      rs = tree.GetNode(treenode->leftChild->rightSib->GetNumber());
      //std::cout << "currndnum: " << currndnum << "\n";
      assert(lc);
      assert(rs);
      assert(treenode->leftChild);
      assert(treenode->leftChild->rightSib);
      //std::cout << "Got here1\n";
      copyNode(treenode->leftChild, lc);
      copyNode(treenode->leftChild->rightSib, rs);

      //std::cout << "Got here2\n";
      nodecopy->leftChild = lc;
      nodecopy->leftChild->rightSib = rs;

      lc->parent = nodecopy;
      rs->parent = nodecopy;

      //std::cout << "Got here3\n";

      // function recursively called on both children
      copyTree(tree, treenode->leftChild, lc);
      copyTree(tree, treenode->leftChild->rightSib, rs);
    } else return;
}

void Tree::copyTopology(const Tree &other){
  //int currNdNum = other.GetNumLeaves();
  Node * r = this->GetNode(other.GetRoot()->GetNumber());
  this->SetRoot(r);
  other.TreeDebug();
  std::cout << "Calling copyTree\n";
  copyTree(*this, other.GetRoot(), r);
}

// assertions fail if full tree is not properly assembled/initialized
void Tree::TreeDebug() const {
    Tree *tr = const_cast<Tree*>(this);
    // need a const node iterator but this will work for now
    Node * r = tr->GetRoot();
    // verify that root's parent points to null
    assert(!r->parent && "r->parent not null");
    int nodeCount = 1;
    assert(r->GetEdgeLen() && "edge length is nonzero");
    PostorderForNodeIterator ptrav = postorder(r);
    Arc arc = ptrav.get();
    while(arc.toNode) {
      nodeCount++;
      //std::cout << "Node encountered: " << arc.fromNode->GetNumber() << "\n";
      assert(arc.toNode->GetEdgeLen() && "edge length is nonzero");
      if (tr->isLeaf2(arc.toNode)) {
          // leaves should not have children
          assert((!arc.toNode->leftChild) && "leaf has no children");
      } else {
          // should have both left and right child
          assert(arc.toNode->leftChild && arc.toNode->leftChild->rightSib && "internal node has two children");
      }
      arc = ptrav.next();
    }
    //std::cout << "Nodes counted: " << nodeCount << "\n";
    //std::cout << "Expected value: " << nodes.size() << "\n";
    assert(nodeCount == nodes.size() && "node count is correct");
}

// look for node with id in tree rooted at n
bool Tree::isNodeConnected(Node *n, int id){
  //std::cout << "Checking whether node " << id << " is connected\n";
  //std::cout << "At " << n->GetNumber() << "\n";

  if(id == n->GetNumber()) return true;
  if(n->IsLeaf()) return false;
  return (this->isNodeConnected(n->leftChild, id) ||
          this->isNodeConnected(n->leftChild->rightSib, id));
  // default case
  std::cout << "Something went wrong in Tree::isNodeConnected(...)\n";
  return false;
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
