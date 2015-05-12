#if !defined(__OPTIMIZE_BRANCHES_H__)
#define __OPTIMIZE_BRANCHES_H__
#include "mt_log.h"
namespace mt {
class MTInstance;
class Arc;
void optimizeAllBranchLengths(MTInstance &cm);
double optimizeAllLengthsForOneEdge(MTInstance & mtInstance, Arc & edge);

} // namespace
#endif
