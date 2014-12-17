#if !defined(__TREE_H__)
#define __TREE_H__
#include <vector>
namespace mt {

typedef unsigned int char_state_t;
class CharStateToPrimitiveInd {
    public:

};
class LeafCharacterVector {
    public:
        LeafCharacterVector()
            :cs2pi(nullptr) {
        }
        LeafCharacterVector(const char_state_t *inp, unsigned len, const CharStateToPrimitiveInd * stateToPrimStates)
            :charVec(inp, inp + len),
            cs2pi(stateToPrimStates) {
        }
        std::vector<char_state_t> charVec;
        const CharStateToPrimitiveInd * cs2pi;
};

typedef std::vector<LeafCharacterVector> CharMatrix;

class PartitionedMatrix {
    public:
        PartitionedMatrix(unsigned numTaxa, const std::vector<unsigned> &numCharsPerPartition)
            :nTaxa(numTaxa),
             nCharsVec(numCharsPerPartition),
             partitions(numCharsPerPartition.size()) {
             for (auto i : partitions) {
                i.resize(numTaxa);
             }
        }
    private:
        unsigned nTaxa;
        std::vector<unsigned> nCharsVec;
        std::vector<CharMatrix> partitions;

};
class ModelDescription {
    public:
        enum AscBiasMode {
            NO_ASC_BIAS = 0,
            VAR_ONLY_NO_MISSING_ASC_BIAS = 1,
            VAR_ONLY_MISSING_ASC_BIAS = 2,
            PARS_ONLY_NO_MISSING_ASC_BIAS = 3,
            PARS_ONLY_MISSING_ASC_BIAS = 4
        };
        ModelDescription(AscBiasMode m)
            :ascBiasMode(m) {
        }
        const AscBiasMode & GetAscBiasMode() const {
            return this->ascBiasMode;
        }
    private:
        AscBiasMode ascBiasMode;
};

}
#endif
