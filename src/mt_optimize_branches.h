#if !defined(__OPTIMIZE_BRANCHES_H__)
#define __OPTIMIZE_BRANCHES_H__
namespace mt {

void optimizeAllBranchLengths(PartitionedMatrix &partMat,
                              Tree &tree,
                              CharModel &cm);


} // namespace
#endif
