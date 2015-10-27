#if ! defined(MT_LIKELHOOD_H)
#define MT_LIKELHOOD_H

namespace mt {
class PartitionedMatrix;
class Tree;
class MTInstance;
class CharModel;

double ScoreTree(PartitionedMatrix &partMat, Tree &tree, MTInstance &instance);
double ScoreTreeForPartition(PartitionedMatrix &partMat, Tree &tree, CharModel &cm, unsigned model);

} // namespace mt

#endif
