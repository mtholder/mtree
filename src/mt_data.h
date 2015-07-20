#if !defined(__MT_DATA_H__)
#define __MT_DATA_H__
#include <vector>
#include "mt_log.h"
namespace mt {

// Some Variables - should these go in a more general header file?
#define MAX_ITERS               100
#define MAX_REARRANGEMENTS      100
#define BZ_EPSILON              1.e-5
#define BRENT_VAR               0.3819660
#define GOLDEN_RAT              1.618034
#define BRAK_VAR1               1.e-20
#define BRAK_VAR2               100.0
#define MT_RATE_MIN             0.0000001
#define MT_RATE_MAX             1000000.0

// Some macros
#define MTREE_SIGN(a,b)         ((b) > 0.0 ? fabs(a) : -fabs(a))

typedef unsigned int char_state_t;
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
            return stateCodeToStateCodeVec.size();
        }
    private:
        std::vector<std::vector<char_state_t> > stateCodeToStateCodeVec;
};

class LeafCharacterVector {
    public:
        LeafCharacterVector()
            :cs2pi(nullptr) {
        }
        LeafCharacterVector(const char_state_t *inp, unsigned len, const CharStateToPrimitiveInd * stateToPrimStates)
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
                          const std::vector<unsigned> &numCharsPerPartition,
                          const std::vector<unsigned> &orig2compressed,
                          const std::vector<double> &patternWts)
            :nTaxa(numTaxa),
             nCharsVec(numCharsPerPartition),
             partitions(numCharsPerPartition.size()),
             origInd2CompressedInd(orig2compressed),
             patternWeights(patternWts) {
             for (auto i : partitions) {
                i.resize(numTaxa);
             }
        }
        unsigned GetNumPartitions() const {
            return this->partitions.size();
        }
        const LeafCharacterVector * GetLeafCharacters(unsigned partIndex, unsigned leafIndex) const {
            return &(this->partitions[partIndex][leafIndex]);
        }
        void fillPartition(unsigned partIndex, const char_state_t ** mat, const CharStateToPrimitiveInd *cs2pi) {
            CharMatrix & part = partitions[partIndex];
            part.resize(nTaxa);
            for (auto i = 0U; i < nTaxa; ++i) {
                part[i] = LeafCharacterVector(mat[i], nCharsVec[partIndex], cs2pi);
            }
        }
    private:
        unsigned nTaxa;
        std::vector<unsigned> nCharsVec;
        std::vector<CharMatrix> partitions;
        std::vector<unsigned> origInd2CompressedInd;
    public:
        const std::vector<double> patternWeights;
};


} //namespace
#endif
