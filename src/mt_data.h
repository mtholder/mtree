#if !defined(__MT_DATA_H__)
#define __MT_DATA_H__
#include <vector>
#include "mt_log.h"
namespace mt {

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
        std::size_t GetNumPartitions() const {
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
