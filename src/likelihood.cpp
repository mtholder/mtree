#include "mt_tree.h"
namespace mt {

void doAnalysis(Tree &tree, CharModel &cm)
{
    NodeIterator *pnit = postorder(tree.GetRoot());
    try {
        Node * c = pnit->get();
        assert(c);
        unsigned partIndex = 0;
        while (c) {
            std::cout << c->number<< "\n";
            const double edgeLen = c->GetEdgeLen();
            if (c->IsLeaf()) {
                const LeafCharacterVector * data = (const LeafCharacterVector *) c->GetData(partIndex);
                LeafWork * work = (LeafWork *) c->GetWork(partIndex);
                cm.fillLeafWork(data, work, edgeLen);
            } else {
                Node * left = c->leftChild;
                Node * n = left->rightSib;
                const double * c1 = left->GetCLA(partIndex);
                const double * c2 = n->GetCLA(partIndex);
                InternalNodeWork * work = (InternalNodeWork *) c->GetWork(partIndex);
                cm.fillInternalWork(c1, c2, work, edgeLen);
            }
            c = pnit->next();
        }
    } catch(...) {
        delete pnit;
        throw;
    }
    delete pnit;

}


void CharModel::fillLeafWork(const LeafCharacterVector *, LeafWork *, double edgeLen) {
}
void CharModel::fillInternalWork(const double * cla1, const double *cla2, InternalNodeWork *, double edgeLen) {
}
}

