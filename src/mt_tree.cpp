#include "mt_tree.h"
#include "mt_optimize_branches.h"

namespace mt {


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
