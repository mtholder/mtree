#if !defined(__MT_INSTANCE_H__)
#define __MT_INSTANCE_H__
#include "mt_tree.h"
#include "mt_data.h"
namespace mt {

class NCL2MT;

enum ProcessActionsEnum {
    SCORE_ACTION=0,
    TREE_SEARCH=1,
    OPTIMIZE_EDGE_LENGTH=2
};

enum TipCaseEnum {
    PLL_TIP_TIP = 0,
    PLL_TIP_INNER = 1,
    PLL_INNER_INNER = 2
};

struct OptimizationSettings {
    unsigned maxIterBrLenSmoothing;
    double brLenDiffThreshold;
    std::vector<bool> partitionConverged;
    std::vector<bool> partitionSmoothed;
    OptimizationSettings()
        :maxIterBrLenSmoothing(20) {
    }
};
bool isTip(int n, std::size_t mxtips);
class TraversalInfo {
    public:
    TipCaseEnum tipCase;
    int pNumber;                  /**< should exist in some nodeptr p->number */
    int qNumber;                  /**< should exist in some nodeptr q->number */
    int rNumber;                  /**< should exist in some nodeptr r->number */
    std::vector<double> qz;
    std::vector<double> rz;
    /* recom */
    int slot_p;                   /**< In recomputation mode, the RAM slot index for likelihood vector of node p, otherwise unused */
    int slot_q;                   /**< In recomputation mode, the RAM slot index for likelihood vector of node q, otherwise unused */
    int slot_r;                   /**< In recomputation mode, the RAM slot index for likelihood vector of node r, otherwise unused */
};

class ConcatModelInfo {
    public:
        PartModelInfo & GetPartitionScoreInfo(std::size_t i) {
            return psiVec.at(i);
        }
        std::vector<PartModelInfo> & GetPartitionScoreInfoVec() {
            return psiVec;
        }
        const PartModelInfo & GetPartitionScoreInfo(std::size_t i) const {
            return psiVec.at(i);
        }
        const std::vector<PartModelInfo> & GetPartitionScoreInfoVec() const {
            return psiVec;
        }
    private:
        std::vector<PartModelInfo > psiVec;
};

class TraversalDescriptor {
    public:
    void SetExecuteMask(std::size_t modelIndex, bool val) {
        executeMask.at(modelIndex) = val;
    }
    void SetParamValue(std::size_t modelIndex, double val) {
        paramV.at(modelIndex) = val;
    }
    void SetQNumber(std::size_t q) {
        qNumber = q;
    }
    std::size_t GetQNumber() const {
        return qNumber;
    }
    void SetPNumber(std::size_t q) {
        pNumber = q;
    }
    std::size_t GetPNumber() const {
        return pNumber;
    }
    TraversalInfo & Top() {
        return topTraversalInfo;
    }
    const TraversalInfo & Top() const {
        return topTraversalInfo;
    }
    private:
    TraversalInfo topTraversalInfo;
    std::deque<bool> executeMask;
    std::vector<double> paramV;
    std::size_t pNumber;
    std::size_t qNumber;
};

class MTInstance {
    friend class NCL2MT;
    ConcatModelInfo si;
    TraversalDescriptor td;
    public:
        PartitionedMatrix partMat;
        Tree tree;
        OptimizationSettings optSettings;
        PartModelInfo & GetPartitionScoreInfo(std::size_t i) {
            return si.GetPartitionScoreInfo(i);
        }
        ConcatModelInfo & GetConcatModelInfo() {
            return si;
        }
        TraversalDescriptor & GetTraversalDescriptor() {
            return td;
        }
        bool HasSearchConverged;
        CharModel & GetCharModel() {
            return *charModelPtr;
        }
        std::size_t GetNumPartitions() const {
            return partMat.GetNumPartitions();
        }
        std::size_t GetNumBranchLenPartitions() const {
            return partScoreInfoByBrLenPart.size();
        }
        std::vector<PartModelInfo *> & GetPartitionScoreInfoForBranchLenPart(std::size_t brLenPartIndex) {
            return partScoreInfoByBrLenPart.at(brLenPartIndex);
        }
        std::vector<PartModelInfo> & GetPartitionScoreInfoVec() {
            return si.GetPartitionScoreInfoVec();
        }
        bool DoAscBiasCorr(std::size_t modelIndex) const;
        std::size_t GetMxTips() const {
            return mxtips;
        }
        bool GetSaveMemory() const {
            return false;
        }
    private:
        MTInstance(const MTInstance &) = delete;
        MTInstance & operator=(const MTInstance &) = delete;
        CharModel * charModelPtr;
        std::vector<std::vector<PartModelInfo*> > partScoreInfoByBrLenPart;
        MTInstance(unsigned numTaxa,
                   const std::vector<unsigned> &numCharsPerPartition,
                   const std::vector<unsigned> &orig2compressed,
                   const std::vector<double> &patternWts,
                   CharModel *cm)
            :partMat(numTaxa, numCharsPerPartition, orig2compressed, patternWts),
            tree(2*numTaxa - 1, numTaxa),
            HasSearchConverged(false),
            charModelPtr(cm) {
        }
        int threadID;
        std::size_t mxtips;
};

void doAnalysis(std::ostream * os,
                MTInstance & mtInstance,
                ProcessActionsEnum action);

void SetExecuteMaskForBranchLengthPart(MTInstance & instance, std::size_t brLenPartIndex, bool val);


inline bool MTInstance::DoAscBiasCorr(std::size_t modelIndex) const {
    const PartModelInfo & psi{si.GetPartitionScoreInfo(modelIndex)};
#   if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
        return psi.CorrectForAscBias() && threadID == 0;
#   else
        return psi.CorrectForAscBias();
#   endif
}

inline void SetExecuteMaskForBranchLengthPart(MTInstance & instance, std::size_t brLenPartIndex, bool val) {
    auto & parts = instance.GetPartitionScoreInfoForBranchLenPart(brLenPartIndex);
    for (auto p : parts) {
        p->SetExecuteMask(val);
    }
}

template<typename T>
inline void SetExecuteMaskByBranchLengthPart(MTInstance & instance, const T & perBranchBool, const bool negate) {
    const auto numBranchParts = instance.GetNumBranchLenPartitions();
    for (auto i = 0U; i < numBranchParts ; ++i) {
        const bool execMaskVal = (negate ? !perBranchBool[i] : perBranchBool[i]);
        SetExecuteMaskForBranchLengthPart(instance, i, execMaskVal);
    }
}

inline void SetExecuteMaskForAllPart(MTInstance & instance, const bool val) {
    auto & psiVec = instance.GetPartitionScoreInfoVec();
    for(auto & psi : psiVec) {
        psi.SetExecuteMask(val);
    }
}

inline void storeExecuteMaskInTraversalDescriptor(TraversalDescriptor & td,
                                                  const PartModelInfoVec & psi) {
    for (auto model = 0U; model < psi.size(); model++) {
        td.SetExecuteMask(model, psi[model].GetExecuteMask());
    }
}

inline void storeValuesInTraversalDescriptor(TraversalDescriptor & td,
                                                  const std::vector<double> & vv) {
    for (auto model = 0U; model < vv.size(); model++) {
        td.SetParamValue(model, vv[model]);
    }
}

inline void storeExecuteMaskInTraversalDescriptor(MTInstance & instance) {
    storeExecuteMaskInTraversalDescriptor(instance.GetTraversalDescriptor(),
                                          instance.GetPartitionScoreInfoVec());
}
inline void storeValuesInTraversalDescriptor(MTInstance & instance,
                                                  const std::vector<double> & vv) {
    storeValuesInTraversalDescriptor(instance.GetTraversalDescriptor(), vv);
}

constexpr double PLL_ZMIN = 1.0E-15;  // max branch prop. to -log(PLL_ZMIN) (= 34)
constexpr double PLL_ZMAX = 1.0 - 1.0E-6; // min branch prop. to 1.0-zmax (= 1.0E-6) 


} // namespace
#endif
