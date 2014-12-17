#include "mt_tree.h"
namespace mt {

void doAnalysis(Tree &tree, CharModel &cm)
{
    NodeIterator *pnit = postorder(tree.GetRoot());
    try {
        Node * c = pnit->get();
        assert(c);
        while (c) {
            std::cout << c->number<< "\n";
            c = pnit->next();
        }
    } catch(...) {
        delete pnit;
        throw;
    }
    delete pnit;

}


}

