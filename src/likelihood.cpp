#include "mt_tree.h"
#include "mt_tree_traversal.h"
#include "mt_char_model.h"
#include <algorithm>
using namespace std;
namespace mt {
void pruneProductStep(const vector<const double *> & v, double * dest, unsigned n) {
    for (auto i = 0U; i < n; ++i) {
        dest[i] = 1.0;
        for (auto c : v) {
            dest[i] *= c[i];
        }
    }
}
void doAnalysis(PartitionedMatrix &partMat, Tree &tree, CharModel &cm)
{
    Node * virtRoot = tree.GetRoot();
    PostorderForNodeIterator pnit = postorder(virtRoot);
    Arc c = pnit.get();
    assert(c.toNode);
    unsigned partIndex = 0;
    unsigned numChars =  c.GetNumChars(partIndex);
    while (c.toNode) {
        std::cout << c.fromNode->number<< "\n";
        const double edgeLen = c.GetEdgeLen();
        if (c.IsFromLeaf()) {
            const LeafCharacterVector * data = c.GetFromNdData(partIndex);
            LeafWork * lw = c.GetFromNdLeafWork(partIndex);
            double * claElements = lw->GetCLAElements();
            double * cla = c.GetFromNdCLA(partIndex, true);
            cm.fillLeafWork(data, claElements, cla, edgeLen, numChars);
        } else {
            vector<const double *> p = c.GetPrevCLAs(partIndex);
            double * beforeArc = c.GetFromNdCLA(partIndex, false);
            pruneProductStep(p, beforeArc, c.GetLenCLA(partIndex));
            double * afterArc = c.GetFromNdCLA(partIndex, true);
            cm.conditionOnSingleEdge(beforeArc, afterArc, edgeLen, numChars);
        }
        c = pnit.next();
    }
    vector<const double *> p = GetSurroundingCLA(virtRoot, nullptr, partIndex);
    double * beforeArc = c.GetFromNdCLA(partIndex, false); // not valid arc, but this should work
    pruneProductStep(p, beforeArc, c.GetLenCLA(partIndex));
    const double lnL = cm.sumLnL(beforeArc, &(partMat.patternWeights[0]), numChars);
    cout << "lnL = " << lnL << "\n";
}
/*
inline void assign(double value, double * dest, unsigned n) {
    for (auto i = 0U; i < n; ++i) {
        dest[i] = value;
    }
}
*/

template<typename T>
inline void assign(T value, T * dest, unsigned n) {
    for (auto i = 0U; i < n; ++i) {
        dest[i] = value;
    }
}


void CharModel::fillLeafWork(const LeafCharacterVector *data,
                             double *claElements,
                             double *cla,
                             double edgeLen,
                             unsigned numChars) {
    /* fill the summed probabilities for each state code */

    const double * tiprob = this->calcTransitionProb(edgeLen);
    const CharStateToPrimitiveInd * s2pi = data->cs2pi;
    const unsigned numStateCodes = s2pi->GetNumStateCodes();
    const unsigned lenCLAWord = nStates*nRateCats;
    double * summedLoc = claElements;
    const unsigned nssq = nStates * nStates;
    assign(0.0, summedLoc, lenCLAWord*numStateCodes);
    for (auto sci = 0U; sci < numStateCodes; ++sci) {
        for (auto toState : s2pi->GetStateCodes(sci)) {
            for (auto ri = 0U; ri < nRateCats; ++ri) {
                for (auto fromState = 0U ; fromState < nStates; ++fromState) {
                    summedLoc[ri*nStates + fromState] += tiprob[ri*nssq + fromState*nStates + toState];
                }
            }
        }
        summedLoc += lenCLAWord;
    }
    /* fill in the cla vector by copying sums */
    for (auto ci = 0U; ci < numChars; ++ci) {
        const char_state_t sc = data->charVec[ci];
        const double * s = claElements + sc * lenCLAWord ;
        copy(s, s + lenCLAWord, cla + ci*lenCLAWord);
    }
}

void CharModel::conditionOnSingleEdge(const double * beforeEdge, double * afterEdge, double edgeLen, unsigned numChars) {
    const double * tiprob = this->calcTransitionProb(edgeLen);
    const unsigned lenCLAWord = nStates*nRateCats;
    const unsigned nssq = nStates * nStates;
    for (auto c = 0U ; c < numChars; ++c) {
        for (auto ri = 0U; ri < nRateCats; ++ri) {
            for (auto fromState = 0U ; fromState < nStates; ++fromState) {
                double prob = 0.0;
                for (auto toState = 0U ; toState < nStates; ++toState) {
                    prob += tiprob[ri*nssq + fromState*nStates + toState]*beforeEdge[ri*nStates + toState];
                }
                afterEdge[ri*nStates + fromState] = prob;
            }
        }
        afterEdge += lenCLAWord;
        beforeEdge += lenCLAWord;
    }
}

double CharModel::sumLnL(const double *cla, const double * patternWeight, unsigned numChars) const {
    const unsigned lenCLAWord = nStates*nRateCats;
    const double * rateCatProb = GetRateCatProb();
    const double * stateFreq = GetRootStateFreq();
    double lnL = 0.0;
    for (auto i = 0U; i < numChars; ++i) {
        double charL = 0.0;
        for (auto ri = 0U; ri < nRateCats; ++ri) {
            double rateL = 0.0;
            for (auto fromState = 0U ; fromState < nStates; ++fromState) {
                rateL += cla[ri*nStates + fromState]*stateFreq[fromState];
            }
            charL += rateL*rateCatProb[ri];
        }
        lnL += patternWeight[i]*log(charL);
        cla += lenCLAWord;
    }
    return lnL;
}


/*
void CharModel::fillInternalWork(const double * cla1, const double *cla2, InternalNodeWork *, double edgeLen) {
//    const double * this->calcTransitionProb(edgeLen);
}
*/


} //namespace
