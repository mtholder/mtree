#include "mt_tree.h"
namespace mt {

void doAnalysis(Tree &tree, CharModel &cm)
{
    std::cout << cm.sumLnL(tree.GetRoot()) << "\n";
}


}

