#include "mt_data.h"
#include <cassert>
namespace mt {

LeafCharacterVector::LeafCharacterVector(
        const char_state_t *inp,
        std::size_t len,
        const CharStateToPrimitiveInd * stateToPrimStates)
    :charVec(len),
    cs2pi(stateToPrimStates) {
    for (auto i = 0U; i < len; ++i) {
        const auto c = inp[i];
        if (getenv("DEBUG_MT_TREE") != nullptr) {
            if (c >= stateToPrimStates->size()) {
                //std::cerr << " " << c << " is too big!\n";
            } else {
                const auto & sc =  stateToPrimStates->GetStateCodes(c);
                //std::cerr << " " << c << " ==> {";
                for (auto s : sc) {
                    //std::cerr << s << ", ";
                }
                //std::cerr << "}\n";
            }
        }
        charVec[i] = inp[i];
    }
    //_DEBUG_VEC(charVec);
}

PartitionedMatrix::PartitionedMatrix(
        unsigned numTaxa,
        const std::vector<double> &patternWts,
        const std::vector<std::size_t> &orig2compressed,
        const std::map<unsigned, std::set<unsigned> > & numStates2PatternIndexSet,
        const std::map<unsigned, std::size_t> & numStates2NumBogusChar)
    :nTaxa(numTaxa),
     nCharsVec(numStates2PatternIndexSet.size()),
     partitions(numStates2PatternIndexSet.size()),
     origInd2CompressedInd(orig2compressed),
     nStatesVec(numStates2PatternIndexSet.size()),
     patternWeights(patternWts) {
    std::size_t partIndex = 0;
    for (auto i : numStates2PatternIndexSet) {
        nCharsVec[partIndex] = i.second.size() + numStates2NumBogusChar.at(i.first);
        nStatesVec[partIndex] = i.first;
        ++partIndex;
    }
    for (auto j : partitions) {
        j.resize(numTaxa);
    }
}

void PartitionedMatrix::fillPartition(unsigned partIndex,
                           const std::vector< std::vector<mt::char_state_t> > & mat,
                           const CharStateToPrimitiveInd *cs2pi) {
    CharMatrix & part = partitions.at(partIndex);
    part.resize(nTaxa);
    if (getenv("DEBUG_MT_TREE") != nullptr) {std::cerr << "PartitionedMatrix::fillPartition \n";}
    for (auto i = 0U; i < nTaxa; ++i) {
        const auto & row = mat[i];
        if (getenv("DEBUG_MT_TREE") != nullptr) {
            //std::cerr << "  taxon " << i << " data:  ";
            for (auto c : row) {
                //std::cerr << c << " ";
            }
            //std::cerr << "\n";
        }
        assert(nCharsVec[partIndex] == row.size());
        part[i] = LeafCharacterVector(&(row[0]), nCharsVec[partIndex], cs2pi);
    }
}

} // namespace mt
