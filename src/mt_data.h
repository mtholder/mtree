#if !defined(__MT_DATA_H__)
#define __MT_DATA_H__
#include <vector>
#include "mt_log.h"
namespace mt {

// Some Variables
#define MAX_ITERS                               100
#define MAX_REARRANGEMENTS                      100
#define BZ_EPSILON                              1.e-5
#define BRENT_VAR                               0.3819660
#define GOLDEN_RAT                              1.618034
#define BRAK_VAR1                               1.e-20
#define BRAK_VAR2                               100.0
#define MT_RATE_MIN                             0.0000001
#define MT_RATE_MAX                             1000000.0
#define MT_ALPHA_MIN                            0.2

// Some macros
#define MTREE_SIGN(a,b)                        ((b) > 0.0 ? fabs(a) : -fabs(a))
#define MT_POINT_GAMMA(prob,alpha,beta)        PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define GetPatData(ind) instance.GetCharModel(ind)

typedef unsigned int char_state_t;
class CharModel;
typedef std::vector<CharModel*> ModelVec;

class CharStateToPrimitiveInd {
    public:
        CharStateToPrimitiveInd(unsigned numStateCodes)
            :stateCodeToStateCodeVec(numStateCodes) {
        }
        const std::vector<char_state_t> & GetStateCodes(unsigned i) const {
            return stateCodeToStateCodeVec[i];
        }
        void SetStateCode(unsigned i, const std::vector<char_state_t> & v) {
            this->stateCodeToStateCodeVec[i] = v;
        }
        unsigned GetNumStateCodes() const {
            return static_cast<unsigned>(stateCodeToStateCodeVec.size());
        }
    private:
        std::vector<std::vector<char_state_t> > stateCodeToStateCodeVec;
};

class LeafCharacterVector {
    public:
        LeafCharacterVector()
            :cs2pi(nullptr) {
        }
        LeafCharacterVector(const char_state_t *inp, std::size_t len, const CharStateToPrimitiveInd * stateToPrimStates)
            :charVec(len),
            cs2pi(stateToPrimStates) {
            for (auto i = 0U; i < len; ++i) {
                charVec[i] = inp[i];
            }
            _DEBUG_VEC(charVec);
        }
        std::vector<char_state_t> charVec;
        const CharStateToPrimitiveInd * cs2pi;
};

typedef std::vector<LeafCharacterVector> CharMatrix;

class PartitionedMatrix {
    public:
        PartitionedMatrix(unsigned numTaxa,
                          const std::vector<std::size_t> &numCharsPerPartition,
                          const std::vector<std::size_t> &orig2compressed,
                          const std::vector<unsigned> &numStatesPerPartition,
                          const std::vector<double> &patternWts)
            :nTaxa(numTaxa),
             nCharsVec(numCharsPerPartition),
             partitions(numCharsPerPartition.size()),
             origInd2CompressedInd(orig2compressed),
             nStatesVec(numStatesPerPartition),
             patternWeights(patternWts) {
             for (auto i : partitions) {
                i.resize(numTaxa);
             }
        }
        std::size_t GetNumPartitions() const {
            return this->partitions.size();
        }
        const LeafCharacterVector * GetLeafCharacters(unsigned partIndex, unsigned leafIndex) const {
            return &(this->partitions[partIndex][leafIndex]);
        }
        void fillPartition(unsigned partIndex, const char_state_t ** mat, const CharStateToPrimitiveInd *cs2pi) {
            CharMatrix & part = partitions.at(partIndex);
            part.resize(nTaxa);
            for (auto i = 0U; i < nTaxa; ++i) {
                part[i] = LeafCharacterVector(mat[i], nCharsVec[partIndex], cs2pi);
            }
        }
    private:
        unsigned nTaxa;
        std::vector<std::size_t> nCharsVec;
        std::vector<CharMatrix> partitions;
        std::vector<std::size_t> origInd2CompressedInd;
        std::vector<unsigned> nStatesVec;
    public:
        const std::vector<double> patternWeights;
};


} //namespace
#endif
