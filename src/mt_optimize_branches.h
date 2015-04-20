#if !defined(__OPTIMIZE_BRANCHES_H__)
#define __OPTIMIZE_BRANCHES_H__
namespace mt {

class MTInstance;
double optimizeAllBranchLengths(MTInstance &cm);
double optimizeAllLengthsForOneEdge(Arc edge, MTInstance & mtInstance);


} // namespace
#endif
