#if !defined(__PAT_CLASS__)
#define __PAT_CLASS__

#include "mt_log.h"
#include "mt_tree.h"
#include "ncl/nxsallocatematrix.h"
#include "ncl/nxsmultiformat.h"
#include <vector>
#include <map>
#include <cassert>
#include <iostream>
#include <cmath>
#include <string>
#include <utility>

namespace mt {

typedef double ** TiMat;
typedef double *** TiMatVec;
typedef void (* TiMatFunc)(double, TiMatVec);
typedef unsigned char BitField;
typedef std::vector<BitField> BitFieldRow;
typedef std::vector<BitFieldRow> BitFieldMatrix;
typedef std::map<BitField, unsigned> BitsToCount;
typedef std::pair<const Node *, unsigned int> NodeID;
typedef std::map<BitField, std::vector<double> > MaskToProbsByState;
typedef std::map<BitField, MaskToProbsByState > MaskToMaskToProbsByState;

const unsigned MAX_NUM_STATES = 8*sizeof(BitField);

class NodeInfo;
class MTInstance;

typedef std::pair<BitField, BitField> MaskPair;
typedef std::vector<MaskPair> VecMaskPair;
typedef std::map<BitField, VecMaskPair> MaskToVecMaskPair;
typedef std::vector<VecMaskPair> VMaskToVecMaskPair;
typedef std::vector<int> stateSetContainer;

int convertBitToIndex(int i);
std::vector<int> subsetsContainingGivenState(int fullSet, int givenState);
std::vector<int> subsetsOfGivenSize(int obsStSet, int numBits);
int getNextCommStSet(const int obsStSet, int i);
std::string convertToBitFieldMatrix(const NxsCharactersBlock & charsBlock, BitFieldMatrix & bfMat);
double pclassCalcTransitionProb(int ancIndex, int i, double edgeLen, MTInstance & instance, unsigned model);
double calcProbOfSubtreeForObsStSetAndComm(NodeInfo * subtreeInfo, int ancIndex, int obsBits, int commonStates, double edgeLen,
                                           MTInstance &instance, unsigned model);
double calcProbOfSubtreeForObsStSetNoRepeated(NodeInfo * subtreeInfo, int ancIndex, int obsBits, double edgeLen, MTInstance &instance, unsigned model);
void cleanVirtualEdgeLens(Node * root);
double calcUninformativePatterns(MTInstance & instance, Node * nd, unsigned charIndex, unsigned model);
double addUninformativePatternProbs(MTInstance & instance);



inline const std::vector<double> * getProbsForStatesMask(const MaskToProbsByState *m, const BitField sc) {
    if (m == 0L)
        return 0L;
    MaskToProbsByState::const_iterator scIt = m->find(sc);
    if (scIt == m->end())
        return 0L;
    return &(scIt->second);
}

inline std::vector<double> * getMutableProbsForStatesMask(MaskToProbsByState *m, const BitField sc) {
    if (m == 0L)
        return 0L;
    MaskToProbsByState::iterator scIt = m->find(sc);
    if (scIt == m->end())
        return 0L;
    return &(scIt->second);
}

class ProbForObsStateSet{ //for each state want to set -1 to 1 and all else to 0 (will be either 1 or anything up to NumStates)
    public:
        ProbForObsStateSet(unsigned int numStates) {
            std::vector<double> initialVal(numStates, 0.0);
            noRepeatedState = initialVal;
            probVec.assign(numStates, initialVal);
        }

        std::vector<double> & getProbForCommState(int commState) {
            if(commState == -1)
                return noRepeatedState;
            return probVec.at(commState);
        }
        void write_pv() const {
          for(std::vector<probvec_t>::const_iterator pIt = probVec.begin(); pIt != probVec.end(); pIt++){
            std::vector<double> vec = *pIt;
            for(std::vector<double>::const_iterator vecIt = vec.begin(); vecIt != vec.end(); vecIt++){
              std::cerr << *vecIt << " ";
            }
            std::cerr << "\n";
          }
        }
        /*void writeDebug(std::ostream & o, const CommonInfo & blob) const {
            o << "ProbForObsStateSet{\n  ";
            o << "-1 ";
            for (unsigned i = 0; i < noRepeatedState.size() ; ++i) {
                o << noRepeatedState[i] << " ";
            }
            for (unsigned i = 0; i < probVec.size() ; ++i) {
                o << "\n  " << i << ' ';
                const probvec_t & pvi = probVec[i];
                for (unsigned i = 0; i < pvi.size() ; ++i) {
                    o << pvi[i] << " ";
                }
            }
            o << "}\n";
        } */

    private:
        typedef std::vector<double> probvec_t;
        std::vector<probvec_t> probVec;
        probvec_t noRepeatedState;
};

} //namespace
#endif
