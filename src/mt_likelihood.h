#if !defined(__LIKELIHOOD_H__)
#define __LIKELIHOOD_H__
#include "mt_log.h"
namespace mt {
class PartitionedMatrix;
class Tree;
class CharModel;
double ScoreTree(const PartitionedMatrix &partMat, const Tree &tree, const CharModel &cm);

} // namespace
#endif
