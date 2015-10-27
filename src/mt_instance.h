#if !defined(__MT_INSTANCE_H__)
#define __MT_INSTANCE_H__

#include "mt_tree.h"
#include "mt_char_model.h"
#include "mt_data.h"

class NxsSimpleTree; // TEMP, adapting should be done in main to remove this

namespace mt {

class NCL2MT;

enum ProcessActionsEnum {
    SCORE_ACTION=0,
    TREE_SEARCH=1,
    OPTIMIZE_PARS=2,
    OPTIMIZE_BR_LEN=3
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
        ~MTInstance();
    private:
        MTInstance(const MTInstance &) = delete;
        MTInstance & operator=(const MTInstance &) = delete;
        std::vector<CharModel*> models;

        MTInstance(const std::map<unsigned, std::vector< std::vector<mt::char_state_t> > > & ns2RawMat,
                   const std::map<unsigned, mt::CharStateToPrimitiveInd> & cs2pi,
                   const std::vector<double> &patternWts,
                   const std::vector<std::size_t> &orig2compressed,
                   const std::map<unsigned, std::set<unsigned> > & numStates2PatternIndexSet,
                   const NxsSimpleTree & nxsTree, //
                   const ModelDescription & md);
};
void doAnalysis(std::ostream * os,
                MTInstance & mtInstance,
                ProcessActionsEnum action);

} // namespace
#endif
