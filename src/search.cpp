#include "search.h"
#include "mt_tree.h"

#include <assert.h>

// Implementation of Tree Topology Search Algorithm ported from PLL (searchAlgo.c and topologies.c)

void performSearch (MTInstance &instance, int steps, Node *p) {
  BestTree *treeList; // Need to allocate space?
  for (int i = 0; (!instance.HasSearchConverged); i++) {
    //Optimize branch lengths of tree here
    treeList[i] = BestTree::BestTree(double ScoreTree(instance.partmat, instance.tree, instance.GetCharModel()),
                               Tree instance.tree);
    searchStep(instance);
  }
  return;
}

void searchStep (MTInstance &instance) {

}
