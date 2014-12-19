#if !defined(__MT_INSTANCE_H__)
#define __MT_INSTANCE_H__
#include "mt_tree.h"
#include "mt_data.h"
namespace mt {


enum ProcessActionsEnum {
    SCORE_ACTION=0
};
struct OptimizationSettings {
    unsigned maxIterBrLenSmoothing;

    OptimizationSettings()
        :maxIterBrLenSmoothing(20) {
    }
};
class MTInstance {
    public:
        MTInstance(unsigned numTaxa,
                   const std::vector<unsigned> &numCharsPerPartition,
                   const std::vector<unsigned> &orig2compressed,
                   const std::vector<double> &patternWts,
                   CharModel *cm)
            :partMat(numTaxa, numCharsPerPartition, orig2compressed, patternWts),
            tree(2*numTaxa - 1, numTaxa),
            charModelPtr(cm) {
        }
        PartitionedMatrix partMat;
        Tree tree;
        OptimizationSettings optSettings;
        CharModel & GetCharModel() {
            return *charModelPtr;
        }
    private:
        MTInstance(const MTInstance &) = delete;
        MTInstance & operator=(const MTInstance &) = delete;
        CharModel * charModelPtr;
        
};
void doAnalysis(std::ostream * os,
                MTInstance & mtInstance,
                ProcessActionsEnum action);

} // namespace
#endif