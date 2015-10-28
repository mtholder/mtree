#include "mt_tree.h"
#include "mt_tree_traversal.h"
#include "mt_char_model.h"
#include "mt_log.h"
#include "mt_likelihood.h"
#include "mt_data.h"
#include "mt_instance.h"
#include <algorithm>
using namespace std;
namespace mt {
void pruneProductStep(const vector<const double *> & v, double * dest, std::size_t n);


void pruneProductStep(const vector<const double *> & v, double * dest, std::size_t n) {
    for (auto i = 0U; i < n; ++i) {
        dest[i] = 1.0;
        for (auto c : v) {
            dest[i] *= c[i];
        }
    }
}

// Calculate likelihood for one partition for a tree
double ScoreTreeForPartition(PartitionedMatrix &partMat, Tree &tree, CharModel &cm, unsigned model) {
  Node * virtRoot = tree.GetRoot();
  virtRoot = virtRoot->leftChild->rightSib;
  //Set up a traversal
  PostorderForNodeIterator pnit = postorder(virtRoot);
  Arc c = pnit.get();
  assert(c.toNode);
  std::size_t numChars =  c.GetNumChars(model);
  while (c.toNode) {
    //_DEBUG_VAL(c.fromNode->number);
    const double edgeLen = c.GetEdgeLen();
    if (c.IsFromLeaf()) {
      const LeafCharacterVector * data = c.GetFromNdData(model);
      LeafWork * lw = c.GetFromNdLeafWork(model);
      double * claElements = lw->GetCLAElements();
      double * cla = c.GetFromNdCLA(model, true);
      cm.fillLeafWork(data, claElements, cla, edgeLen, numChars);

      //_DEBUG_CLA(cla, cm.GetNumRates(), cm.GetNumStates(), numChars);
    } else {
      vector<const double *> p = c.GetPrevCLAs(model);
      double * beforeArc = c.GetFromNdCLA(model, false);
      pruneProductStep(p, beforeArc, c.GetLenCLA(model));
      double * afterArc = c.GetFromNdCLA(model, true);
      cm.conditionOnSingleEdge(beforeArc, afterArc, edgeLen, numChars);

      //_DEBUG_CLA(beforeArc, cm.GetNumRates(), cm.GetNumStates(), numChars);
      //_DEBUG_CLA(afterArc, cm.GetNumRates(), cm.GetNumStates(), numChars);
    }
      c = pnit.next();
  }
  // create the cla for the virtual root
  vector<const double *> p = GetSurroundingCLA(virtRoot, nullptr, model);
  double * beforeArc = c.GetFromNdCLA(model, false); // not valid arc, but this should work
  pruneProductStep(p, beforeArc, c.GetLenCLA(model));

  //_DEBUG_CLA(beforeArc, cm.GetNumRates(), cm.GetNumStates(), numChars);
  //_DEBUG_VEC(partMat.patternWeights);

  return cm.sumLnL(beforeArc, &(partMat.patternWeights[0]), numChars);
}

// Calculates Likelihood score for given tree and for all partitions
double ScoreTree(PartitionedMatrix &partMat,
                 Tree &tree,
                 MTInstance &instance, 
                 bool forceRecalc) {
  double result = 0.0;
  //_DEBUG_VAL(result);
   unsigned numParts = static_cast<unsigned>(partMat.GetNumPartitions());
  for(unsigned partIndex = 0; partIndex < numParts; partIndex++){
    //_DEBUG_VAL(instance.dirtyFlags[partIndex]);
    if(forceRecalc || instance.dirtyFlags[partIndex]) {
      instance.likelihoods[partIndex] = ScoreTreeForPartition(partMat,tree,instance.GetCharModel(partIndex),partIndex);
    }
    result += instance.likelihoods[partIndex];
    _DEBUG_FVAL(partIndex); _DEBUG_MVAL(instance.likelihoods[partIndex]); _DEBUG_LVAL(result);
  }
  return result;
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
                             std::size_t numChars) {
    /* fill the summed probabilities for each state code */

    const double * tiprob = this->calcTransitionProb(edgeLen);
    const CharStateToPrimitiveInd * s2pi = data->GetCharStateToPrimitiveInd();
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
    //_DEBUG_CLA(claElements, nRateCats, nStates, numStateCodes);
    /* fill in the cla vector by copying sums */
    for (std::size_t ci = 0U; ci < numChars; ++ci) {
        const char_state_t sc = data->GetCharVec()[ci];
        const double * s = claElements + sc * lenCLAWord ;
        copy(s, s + lenCLAWord, cla + ci*lenCLAWord);
    }
}

void CharModel::conditionOnSingleEdge(const double * beforeEdge, double * afterEdge, double edgeLen, std::size_t numChars) {
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
                         std::size_t numChars) const {
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
            //_DEBUG_VAL(charL);
            //_DEBUG_VAL(log(charL));
        }
        lnL += patternWeight[i]*log(charL);
        //_DEBUG_VAL(lnL);
        cla += lenCLAWord;
    }
    return lnL;
}

double MkVarNoMissingAscCharModel::sumLnL(const double *cla,
                         const double * patternWeight,
                         std::size_t numChars) const {
    std::size_t numRealPatterns =  numChars - nStates;
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
                         std::size_t numChars) const {
    std::size_t numRealPatterns =  numChars - nStates;
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
                         std::size_t numChars) const {
    std::size_t numRealPatterns =  numChars - nStates;
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
                         std::size_t numChars) const {
    std::size_t numRealPatterns =  numChars - nStates;
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
