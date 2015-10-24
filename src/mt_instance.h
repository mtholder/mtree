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
        bool curvatOK;
        CharModel & GetCharModel(int i) const {
          return *(models[i]);
        }
        std::vector<CharModel*> GetModelVec() const {
          return models;
        }
        /*void changeRate(int pos, double val) {
          this->charModelPtr->alterRateFreq(pos, val);
        }*/
        double curLikelihood;
        unsigned numPartitions;
        std::vector<double> likelihoods; // stores likelihoods for partitions
        std::vector<bool> dirtyFlags;    // vector of bools to indicate whether partition likelihoods have recently changed (true if dirty)

    private:
        MTInstance(const MTInstance &) = delete;
        MTInstance & operator=(const MTInstance &) = delete;
        std::vector<CharModel*> models;

        MTInstance(unsigned numTaxa,
                   const std::vector<std::size_t> &numCharsPerPartition,
                   const std::vector<unsigned> &numStatesPerPartition,
                   const std::vector<std::size_t> &orig2compressed,
                   const std::vector<double> &patternWts,
                   std::vector<CharModel*> mods)
            :partMat(numTaxa, numCharsPerPartition, orig2compressed, numStatesPerPartition, patternWts),
            tree(2*numTaxa - 1, numTaxa),
            HasSearchConverged(false),
            curvatOK(true),
            numPartitions(static_cast<unsigned>(numCharsPerPartition.size())),
            likelihoods(numCharsPerPartition.size(),0.0),
            dirtyFlags(numCharsPerPartition.size(),false),
            models(mods) {
        }
};
void doAnalysis(std::ostream * os,
                MTInstance & mtInstance,
                ProcessActionsEnum action);

} // namespace
#endif
