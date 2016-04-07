#if ! defined(MT_LIKELHOOD_H)
#define MT_LIKELHOOD_H

namespace mt {
class PartitionedMatrix;
class Tree;
class MTInstance;
class CharModel;

void setScoreFlags(Tree &tree, Node *nd);
void resetScoreFlags(Tree &tree);
double ScoreTree(PartitionedMatrix &partMat, Tree &tree, MTInstance &instance, bool forceRecalc);
double ScoreTreeForPartition(PartitionedMatrix &partMat, Tree &tree, CharModel &cm, unsigned model);
double ScoreTreeOneArcDown(PartitionedMatrix &partMat, Tree &tree, MTInstance &instance, Arc dirtyArc);
double ScoreTreeForPartitionOneArcDown(PartitionedMatrix &partMat, Tree &tree, CharModel &cm, unsigned model, Arc dirtyArc);

} // namespace mt

#endif
