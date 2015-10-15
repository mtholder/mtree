#if !defined(__MT_INSTANCE_H__)
#define __MT_INSTANCE_H__

#include "mt_tree.h"
#include "mt_char_model.h"
#include "mt_data.h"

namespace mt {

class NCL2MT;

enum ProcessActionsEnum {
    SCORE_ACTION=0,
    TREE_SEARCH=1
};
struct OptimizationSettings {
    unsigned maxIterBrLenSmoothing;
    double brLenDiffThreshold;
    bool partitionConverged[50];
    bool partitionSmoothed[50];
    OptimizationSettings()
        :maxIterBrLenSmoothing(20) {
    }
};


class MTInstance {
    friend class NCL2MT;
    public:
        PartitionedMatrix partMat;
        Tree tree;
        OptimizationSettings optSettings;
        bool HasSearchConverged;
        CharModel & GetCharModel() const {
            return *charModelPtr;
        }
        void changeRate(int pos, double val) {
          this->charModelPtr->alterRateFreq(pos, val);
        }
        double curLikelihood;

    private:
        MTInstance(const MTInstance &) = delete;
        MTInstance & operator=(const MTInstance &) = delete;
        CharModel * charModelPtr;

        MTInstance(unsigned numTaxa,
                   const std::vector<unsigned> &numCharsPerPartition,
                   const std::vector<unsigned> &numStatesPerPartition,
                   const std::vector<unsigned> &orig2compressed,
                   const std::vector<double> &patternWts,
                   CharModel *cm)
            :partMat(numTaxa, numCharsPerPartition, orig2compressed, numStatesPerPartition, patternWts),
            tree(2*numTaxa - 1, numTaxa),
            HasSearchConverged(false),
            charModelPtr(cm) {
        }
};
void doAnalysis(std::ostream * os,
                MTInstance & mtInstance,
                ProcessActionsEnum action);

} // namespace
#endif
