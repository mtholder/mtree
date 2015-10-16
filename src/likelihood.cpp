#include "mt_tree.h"
#include "mt_tree_traversal.h"
#include "mt_char_model.h"
#include "mt_log.h"
#include "mt_data.h"
#include "search.h"
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

double ScoreTree(PartitionedMatrix &partMat, Tree &tree, CharModel &cm) {
    Node * virtRoot = tree.GetRoot();
    virtRoot = virtRoot->leftChild->rightSib;
    //Set up a traversal
    PostorderForNodeIterator pnit = postorder(virtRoot);
    Arc c = pnit.get();
    assert(c.toNode);
    unsigned partIndex = 0;
    unsigned numChars =  c.GetNumChars(partIndex);
    while (c.toNode) {
        _DEBUG_VAL(c.fromNode->number);
        const double edgeLen = c.GetEdgeLen();
        if (c.IsFromLeaf()) {
            const LeafCharacterVector * data = c.GetFromNdData(partIndex);
            LeafWork * lw = c.GetFromNdLeafWork(partIndex);
            double * claElements = lw->GetCLAElements();
            double * cla = c.GetFromNdCLA(partIndex, true);
            cm.fillLeafWork(data, claElements, cla, edgeLen, numChars);

            _DEBUG_CLA(cla, cm.GetNumRates(), cm.GetNumStates(), numChars);
        } else {
            vector<const double *> p = c.GetPrevCLAs(partIndex);
            double * beforeArc = c.GetFromNdCLA(partIndex, false);
            pruneProductStep(p, beforeArc, c.GetLenCLA(partIndex));
            double * afterArc = c.GetFromNdCLA(partIndex, true);
            cm.conditionOnSingleEdge(beforeArc, afterArc, edgeLen, numChars);

            _DEBUG_CLA(beforeArc, cm.GetNumRates(), cm.GetNumStates(), numChars);
            _DEBUG_CLA(afterArc, cm.GetNumRates(), cm.GetNumStates(), numChars);
        }
        c = pnit.next();
    }
    // create the cla for the virtual root
    vector<const double *> p = GetSurroundingCLA(virtRoot, nullptr, partIndex);
    double * beforeArc = c.GetFromNdCLA(partIndex, false); // not valid arc, but this should work
    pruneProductStep(p, beforeArc, c.GetLenCLA(partIndex));

    _DEBUG_CLA(beforeArc, cm.GetNumRates(), cm.GetNumStates(), numChars);
    _DEBUG_VEC(partMat.patternWeights);

    // accumulate the cla into a lnL
    return cm.sumLnL(beforeArc, &(partMat.patternWeights[0]), numChars);
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
    _DEBUG_CLA(claElements, nRateCats, nStates, numStateCodes);
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

double CharModel::sumLnL(const double *cla,
                         const double * patternWeight,
                         unsigned numChars) const {
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
            _DEBUG_VAL(charL);
            _DEBUG_VAL(log(charL));
        }
        lnL += patternWeight[i]*log(charL);
        _DEBUG_VAL(lnL);
        cla += lenCLAWord;
    }
    return lnL;
}

double MkVarNoMissingAscCharModel::sumLnL(const double *cla,
                         const double * patternWeight,
                         unsigned numChars) const {
    unsigned numRealPatterns =  numChars - nStates;
    double uncorrLnL = CharModel::sumLnL(cla, patternWeight, numRealPatterns);
    const double fake = 1.0;
    double oneStateCorrectionLnL = CharModel::sumLnL(cla + numChars + 1 - nStates, &fake, 1);
    double oneStateCorrectionL = exp(oneStateCorrectionLnL);
    double correctionL = 1 - (nStates * oneStateCorrectionL);
    double corrLnL = log(correctionL);
    double sw = 0.0;
    for (auto i = 0U; i < numRealPatterns; ++i) {
        sw = patternWeight[i];
    }
    double totalCorrection = corrLnL*sw;
    return uncorrLnL - totalCorrection;
}

// Placeholder - same as no missing right now
double MkVarMissingAscCharModel::sumLnL(const double *cla,
                         const double * patternWeight,
                         unsigned numChars) const {
    unsigned numRealPatterns =  numChars - nStates;
    double uncorrLnL = CharModel::sumLnL(cla, patternWeight, numRealPatterns);
    const double fake = 1.0;
    double oneStateCorrectionLnL = CharModel::sumLnL(cla + numChars + 1 - nStates, &fake, 1);
    double oneStateCorrectionL = exp(oneStateCorrectionLnL);
    double correctionL = 1 - (nStates * oneStateCorrectionL);
    double corrLnL = log(correctionL);
    double sw = 0.0;
    for (auto i = 0U; i < numRealPatterns; ++i) {
        sw = patternWeight[i];
    }
    double totalCorrection = corrLnL*sw;
    return uncorrLnL - totalCorrection;
}

// Placeholder
double MkParsInfNoMissingModel::sumLnL(const double *cla,
                         const double * patternWeight,
                         unsigned numChars) const {
    unsigned numRealPatterns =  numChars - nStates;
    double uncorrLnL = CharModel::sumLnL(cla, patternWeight, numRealPatterns);
    const double fake = 1.0;
    double oneStateCorrectionLnL = CharModel::sumLnL(cla + numChars + 1 - nStates, &fake, 1);
    double oneStateCorrectionL = exp(oneStateCorrectionLnL);
    double correctionL = 1 - (nStates * oneStateCorrectionL);
    double corrLnL = log(correctionL);
    double sw = 0.0;
    for (auto i = 0U; i < numRealPatterns; ++i) {
        sw = patternWeight[i];
    }
    double totalCorrection = corrLnL*sw;
    return uncorrLnL - totalCorrection;
}

// Placeholder
double MkParsInfMissingModel::sumLnL(const double *cla,
                         const double * patternWeight,
                         unsigned numChars) const {
    unsigned numRealPatterns =  numChars - nStates;
    double uncorrLnL = CharModel::sumLnL(cla, patternWeight, numRealPatterns);
    const double fake = 1.0;
    double oneStateCorrectionLnL = CharModel::sumLnL(cla + numChars + 1 - nStates, &fake, 1);
    double oneStateCorrectionL = exp(oneStateCorrectionLnL);
    double correctionL = 1 - (nStates * oneStateCorrectionL);
    double corrLnL = log(correctionL);
    double sw = 0.0;
    for (auto i = 0U; i < numRealPatterns; ++i) {
        sw = patternWeight[i];
    }
    double totalCorrection = corrLnL*sw;
    return uncorrLnL - totalCorrection;
}

} // namespace

#include "mt_instance.h"

namespace mt {


void doAnalysis(ostream * os, MTInstance & instance, enum ProcessActionsEnum action) {
    action = SCORE_ACTION;
    if (action == SCORE_ACTION) {
        const double lnL = ScoreTree(instance.partMat, instance.tree, instance.GetCharModel());
        if (os) {
            *os << "lnL = " << lnL << "\n";
        }
    } else if (action == TREE_SEARCH) {
      //int steps = 10;
      double startL = ScoreTree(instance.partMat, instance.tree, instance.GetCharModel());
      *os << "Starting likelihood = " << startL << "\n";
      Node * p = instance.tree.GetLeaf(4)->parent->parent;
      mtreeTestSPR (instance, p, 2, startL);
      double endL = ScoreTree(instance.partMat, instance.tree, instance.GetCharModel());
      *os << "Likelihood after subtree removed = " << endL << "\n";
    //  performSearch(instance, steps, instance.tree);
    }
}


} //namespace
