#include "mt_instance.h"
#include "mt_data.h"
#include "mt_tree.h"
#include "mt_tree_traversal.h"
#include "pattern_class.h"
#include "parsimony.h"
#include "ncl/nxsallocatematrix.h"
#include "ncl/nxsmultiformat.h"

#include <utility>
#include <vector>
#include <cassert>

// Functions for calculating pattern class probabilities as implemented in PhyPatClassProb
// by Mark Holder and Jordan Koch
// Used to correct for ascertainment bias and parsimony-informative data
// taken from uninformative_case.cpp

namespace mt {


//void initSettings(MTInstance &instance) {

// free nodeIDToProbInfo data structure for each node
 void freeProbInfo(PostorderForNodeIterator iter, NodeIDToProbInfo & nodeIDToProbInfo) {
      Arc travArc = iter.get();
      while (travArc.toNode)
      {
        const Node * nd = travArc.fromNode;
        std::vector<Node *> children = nd->GetChildren();
        const unsigned numChildren = children.size();
        NodeID currNdID(nd, 0);
        if (numChildren > 0)
          delete nodeIDToProbInfo[currNdID];
        travArc = iter.next();
      }
}

void ProbInfo::createForTip(const MTInstance & instance) {
    this->byParsScore.resize(1);
    ProbForParsScore & forZeroSteps = this->byParsScore[0];
    unsigned stateIndex = 0;
    for (std::vector<BitField>::const_iterator scIt = GetPatData(0).singleStateCodes.begin();
            scIt != GetPatData(0).singleStateCodes.end();
            ++scIt, ++stateIndex) {
        const BitField sc = *scIt;
        std::vector<double> & pVec = forZeroSteps.byDownPass[sc][sc];
        pVec.assign(GetPatData(0).GetNumRates()*GetPatData(0).GetNumStates(), 0.0);
        for (unsigned r = 0; r < GetPatData(0).GetNumRates(); ++r)
            pVec[GetPatData(0).GetNumStates()*r + stateIndex] = 1.0;
    }
    this->nLeavesBelow = 1;
}

/*
void ProbInfo::addToAncProbVecSymmetric(std::vector<double> & pVec,
        const double *** leftPMatVec, const std::vector<double> * leftProbs,
        const double *** rightPMatVec, const std::vector<double> * rightProbs,
        const std::vector<unsigned int> & rightChildStateCodeTranslation,
        const MTInstance & instance) {
    if (leftProbs == 0L || rightProbs == 0L)
        return;
    unsigned rOffset = 0;
    // ignore loop over rates blob.nRates = 1
    for (unsigned r = 0; r < GetPatData(0).GetNumRates(); ++r) {
        // this code looks up the correct transition prob matrix
        const double ** leftPMat = leftPMatVec[r];
        const double ** rightPMat = rightPMatVec[r];


        for (unsigned ancState = 0; ancState < GetPatData(0).GetNumStates(); ++ancState) {
            double leftProb = 0.0;
            double rightProb = 0.0;
            for (unsigned desState = 0; desState < GetPatData(0).GetNumStates(); ++desState) {

                const double leftTiProb = leftPMat[ancState][desState];
                const double leftAccumProb = (*leftProbs)[rOffset + desState];
                leftProb += leftTiProb*leftAccumProb;


                const unsigned int translatedState = rightChildStateCodeTranslation[desState];
                const double rightTiProb = rightPMat[ancState][translatedState];
                const double rightAccumProb = (*rightProbs)[rOffset + desState];
                std::cerr << "addToAncProbVecSymmetric tr = " << translatedState << " des " << desState << " rightTiProb= " << rightTiProb << " rightAccumProb = " << rightAccumProb <<'\n';
                rightProb += rightTiProb*rightAccumProb;
            }
            pVec[rOffset + ancState] += leftProb*rightProb;
        }
        rOffset += GetPatData(0).GetNumStates();
    }
}
*/
/*
void ProbInfo::addToAncProbVec(std::vector<double> & pVec,
        const double *** leftPMatVec, const std::vector<double> * leftProbs,
        const double *** rightPMatVec, const std::vector<double> * rightProbs,
        const MTInstance & instance) {
    if (leftProbs == 0L || rightProbs == 0L)
        return;
    unsigned rOffset = 0;
    // ignore loop over rates blob.nRates = 1
    for (unsigned r = 0; r < GetPatData(0).GetNumRates(); ++r) {
        // this code looks up the correct transition prob matrix
        const double ** leftPMat = leftPMatVec[r];
        const double ** rightPMat = rightPMatVec[r];


        for (unsigned ancState = 0; ancState < GetPatData(0).GetNumStates(); ++ancState) {
            double leftProb = 0.0;
            double rightProb = 0.0;
            for (unsigned desState = 0; desState < GetPatData(0).GetNumStates(); ++desState) {
                const double leftTiProb = leftPMat[ancState][desState];
                const double leftAccumProb = (*leftProbs)[rOffset + desState];
                leftProb += leftTiProb*leftAccumProb;



                const double rightTiProb = rightPMat[ancState][desState];
                const double rightAccumProb = (*rightProbs)[rOffset + desState];
                rightProb += rightTiProb*rightAccumProb;
            }
            pVec[rOffset + ancState] += leftProb*rightProb;
        }
        rOffset += GetPatData(0).GetNumStates();
    }
}

std::set<BitField> toElements(BitField sc) {
    std::set<BitField> ret;
    BitField curr = 1;
    while (curr <= sc) {
        if (curr & sc)
            ret.insert(curr);
        curr *= 2;
    }
    return ret;
}
*/
std::string convertToBitFieldMatrix(const NxsCharactersBlock & charsBlock,
                    BitFieldMatrix & bfMat) {
    std::vector<const NxsDiscreteDatatypeMapper *> mappers = charsBlock.GetAllDatatypeMappers();
    assert(!(mappers.empty() || mappers[0] == NULL));
    // Fix once algorithm allows for missing data/mixed datasets
    assert(mappers.size() == 1);

    NxsUnsignedSet scratchSet;
    NxsUnsignedSet * toInclude;
    for (unsigned i = 0; i < charsBlock.GetNChar(); ++i)
        scratchSet.insert(i);
    toInclude = & scratchSet;
    std::set <const NxsDiscreteDatatypeMapper *> usedMappers;
    for (NxsUnsignedSet::const_iterator indIt = toInclude->begin(); indIt != toInclude->end(); ++indIt) {
        unsigned charIndex = *indIt;
        usedMappers.insert(charsBlock.GetDatatypeMapperForChar(charIndex));
    }
    assert(usedMappers.size() == 1);
    const NxsDiscreteDatatypeMapper & mapper = **usedMappers.begin();
    NxsCharactersBlock::DataTypesEnum inDatatype = mapper.GetDatatype();
    const unsigned nStates =  mapper.GetNumStates();
    //assert(nStates <= MAX_NUM_STATES);
    for (NxsDiscreteStateCell i = 0; i < (NxsDiscreteStateCell)nStates; ++i)
        assert(mapper.GetStateSetForCode(i).size() == 1);
    const std::string fundamentalSymbols = mapper.GetSymbols();
    assert(fundamentalSymbols.length() == nStates);
    const unsigned nTaxa = charsBlock.GetNTax();
    const unsigned includedNChar = toInclude->size();
    bfMat.resize(nTaxa);
    for (unsigned i = 0; i < nTaxa; ++i) {
        BitFieldRow & bfRow = bfMat[i];
        bfRow.resize(includedNChar);
        const NxsDiscreteStateRow & row = charsBlock.GetDiscreteMatrixRow(i);
        assert(!row.empty());
        unsigned j = 0;
        for (NxsUnsignedSet::const_iterator tIncIt = toInclude->begin(); tIncIt != toInclude->end(); ++tIncIt, ++j) {
            const NxsDiscreteStateCell & cell = row.at(*tIncIt);
            // Change assert when accounted for missing data
            assert(!(cell < 0 || cell >= (NxsDiscreteStateCell) nStates));
            const int bfi = 1 << (int) cell;
            bfRow[j] = BitField(bfi);
        }
    }
    return fundamentalSymbols;
}

/*
void initInfo(MTInstance &instance) {
    GetPatData(0).isMkvSymm = false;
    GetPatData(0).pVecLen = GetPatData(0).GetNumStates()*GetPatData(0).GetNumRates();
    GetPatData(0).categStateProb.assign(GetPatData(0).pVecLen, 1.0/((double)GetPatData(0).pVecLen));
    GetPatData(0).singleStateCodes.clear();
    GetPatData(0).multiStateCodes.clear();
    GetPatData(0).stateCodesToSymbols.clear();
    unsigned lbfU = (1 << GetPatData(0).GetNumStates()) - 1;
    GetPatData(0).stateCodeToNumStates.assign(lbfU + 1, 0);
    GetPatData(0).lastBitField = BitField(lbfU);
    BitField sc = 1;
    for (;; sc++) {
        const std::set<BitField> sbf = toElements(sc);
        if (sbf.size() == 1) {
            unsigned stInd = GetPatData(0).singleStateCodes.size();
            GetPatData(0).singleStateCodes.push_back(sc);
            GetPatData(0).stateIndexToStateCode.push_back(sc);
            assert(GetPatData(0).stateIndexToStateCode[stInd] == sc);
            GetPatData(0).stateCodesToSymbols[sc] = GetPatData(0).alphabet[stInd]; // FIX THIS
            GetPatData(0).stateCodeToNumStates.at(sc) = 1;
        } else {
            GetPatData(0).multiStateCodes.push_back(sc);
            std::string sym;
            for (std::set<BitField>::const_iterator sbfIt = sbf.begin(); sbfIt != sbf.end(); ++sbfIt)
                sym.append(GetPatData(0).stateCodesToSymbols[*sbfIt]);
            GetPatData(0).stateCodeToNumStates.at(sc) = sbf.size();
            GetPatData(0).stateCodesToSymbols[sc] = sym;
        }

        if (sc == GetPatData(0).lastBitField)
            break;

    }
    GetPatData(0).pairsForUnionForEachDownPass.clear();
    GetPatData(0).pairsForUnionForEachDownPass.resize(GetPatData(0).lastBitField + 1);
    GetPatData(0).pairsForIntersectionForEachDownPass.clear();
    GetPatData(0).pairsForIntersectionForEachDownPass.resize(GetPatData(0).lastBitField + 1);

    for (sc = 1;; sc++) {
        if (GetPatData(0).stateCodeToNumStates[sc] > 1) {
            VecMaskPair & forUnions = GetPatData(0).pairsForUnionForEachDownPass[sc];
            for (BitField leftSC = 1; leftSC < sc; ++leftSC) {
                if ((leftSC | sc) != sc)
                    continue;
                BitField rightSC = sc - leftSC;
                assert((rightSC | leftSC) == sc);
                forUnions.push_back(MaskPair(leftSC, rightSC));
            }
        }
        VecMaskPair & forIntersections = GetPatData(0).pairsForIntersectionForEachDownPass[sc];
        for (BitField leftSC = 1; leftSC <= GetPatData(0).lastBitField ; ++leftSC) {
            for (BitField rightSC = 1; rightSC <= GetPatData(0).lastBitField ; ++rightSC) {
                if ((leftSC & rightSC) != sc)
                    continue;
                forIntersections.push_back(MaskPair(leftSC, rightSC));
            }
        }

        if (sc == GetPatData(0).lastBitField)
            break;

    }
        GetPatData(0).statesSupersets.clear();
    GetPatData(0).statesSupersets.resize(GetPatData(0).lastBitField + 1);
    for (sc = 1;;++sc) {
        BitFieldRow & ssRow = GetPatData(0).statesSupersets[sc];
        for (BitField ss = sc;; ++ss) {
            if ((ss & sc) == sc)
                ssRow.push_back(ss);
            if (ss == GetPatData(0).lastBitField)
                break;
        }

        if (sc == GetPatData(0).lastBitField)
            break;

    }

    for (sc = 1;; sc++) {
        const VecMaskPair & forUnions = GetPatData(0).pairsForUnionForEachDownPass[sc];
        const VecMaskPair & forIntersections = GetPatData(0).pairsForIntersectionForEachDownPass[sc];
        if (sc == GetPatData(0).lastBitField)
            break;
    }
}
*/

/*
// \returns true if there were probabalities that were summed
bool ProbInfo::allCalcsForAllPairs(
            MaskToProbsByState & forCurrScoreDownPass,
            const VecMaskPair & pairVec,
            const ProbInfo & leftPI,
            const double *** leftPMatVec,
            const ProbInfo & rightPI,
            const double *** rightPMatVec,
            const unsigned accumScore,
            const bool doingIntersection,
            const MTInstance & instance)
{   // order (2^k)^2
    const unsigned leftMaxP = leftPI.getMaxParsScore();
    const unsigned rightMaxP = rightPI.getMaxParsScore();
    bool probsAdded = false;
    for (VecMaskPair::const_iterator fuIt = pairVec.begin(); fuIt != pairVec.end(); ++fuIt) {
        const BitField leftDown = fuIt->first;
        const BitField rightDown = fuIt->second;
        std::cerr << "from line: " << __LINE__<< ": "; std::cerr << "leftDown = " << GetPatData(0).toSymbol(leftDown) << " rightDown = " << GetPatData(0).toSymbol(rightDown) << '\n';
        assert(leftDown > 0);
        assert(rightDown > 0);
        if (doingIntersection) {
            assert((leftDown & rightDown) != 0);
        }
        else {
            assert((leftDown & rightDown) == 0);
        }
        const unsigned leftMinAccum = GetPatData(0).getnstates(leftDown) - 1;
        const unsigned rightMinAccum = GetPatData(0).getnstates(rightDown) - 1;
        if (leftMinAccum + rightMinAccum > accumScore) {
            std::cerr << "from line: " << __LINE__<< ": minScore exceed required score\n";
            continue;
        }
        const unsigned acMaxLeftAccum = std::min(accumScore - rightMinAccum, leftMaxP);
        const unsigned acMaxRightAccum = std::min(accumScore - leftMinAccum, rightMaxP);
        const unsigned acMinLeftAccum = std::max(leftMinAccum, accumScore - acMaxRightAccum);
        if (false) {
                std::cerr << "accumScore = " << accumScore << '\n';
                std::cerr << "accumScore = " << accumScore << '\n';
                std::cerr << "leftMinAccum = " << leftMinAccum << '\n';
                std::cerr << "rightMinAccum = " << rightMinAccum << '\n';
                std::cerr << "acMaxLeftAccum = " << acMaxLeftAccum << '\n';
                std::cerr << "acMaxRightAccum = " << acMaxRightAccum << '\n';
                std::cerr << "acMinLeftAccum = " << acMinLeftAccum << '\n';
        }
        // order N
        for (unsigned leftAccum = acMinLeftAccum; leftAccum <= acMaxLeftAccum; ++leftAccum) {
            const unsigned rightAccum = accumScore - leftAccum;
            assert(rightAccum <= rightMaxP);
            const ProbForParsScore & leftFPS = leftPI.getByParsScore(leftAccum);
            const MaskToProbsByState * leftM2PBS = leftFPS.getMapPtrForDownPass(leftDown);
            if (leftM2PBS == 0L) {
                std::cerr << "from line: " << __LINE__<< ": left child empty row. Skipping...\n";
                continue;
            }
            const ProbForParsScore & rightFPS = rightPI.getByParsScore(rightAccum);
            const MaskToProbsByState * rightM2PBS = rightFPS.getMapPtrForDownPass(rightDown);
            if (rightM2PBS == 0L) {
                std::cerr << "from line: " << __LINE__<< ": right child empty row. Skipping...\n";
                continue;
            }
            const BitFieldRow & leftSSRow = GetPatData(0).statesSupersets[leftDown];
            // order (2^k)
            for (BitFieldRow::const_iterator lasIt = leftSSRow.begin(); lasIt != leftSSRow.end(); ++lasIt) {
                const BitField leftAllStates = *lasIt;
                    std::cerr << "from line: " << __LINE__<< ":  leftAllStates="  << GetPatData(0).toSymbol(leftAllStates) << '\n';
                const std::vector<double> * leftProbs = getProbsForStatesMask(leftM2PBS, leftAllStates);
                if (leftProbs == 0L) {
                    std::cerr << "from line: " << __LINE__<< ": left child empty bin. Skipping...\n";
                    continue;
                }
                const BitFieldRow & rightSSRow = GetPatData(0).statesSupersets[rightDown];
                // order (2^k)
                for (BitFieldRow::const_iterator rasIt = rightSSRow.begin(); rasIt != rightSSRow.end(); ++rasIt) {
                    const BitField rightAllStates = *rasIt;
//                  std::cerr << "from line: " << __LINE__<< ":  rightAllStates="  << blob.toSymbol(rightAllStates) << '\n';
                    const std::vector<double> * rightProbs = getProbsForStatesMask(rightM2PBS, rightAllStates);
                    if (rightProbs == 0L) {
//                      std::cerr << "from line: " << __LINE__<< ": right child empty bin. Skipping...\n";
                        continue;
                    }
                    const BitField ancAllField = (rightAllStates|leftAllStates);

                    std::vector<double> * ancVec = getMutableProbsForStatesMask(&forCurrScoreDownPass, ancAllField);
                    if (ancVec == 0L) {
                        ancVec = &(forCurrScoreDownPass[ancAllField]);
                        ancVec->assign(GetPatData(0).GetNumRates()*GetPatData(0).GetNumStates(), 0.0);
                    }
                    probsAdded = true;
//                  std::cerr << __LINE__ << " adding:";
//                  std::cerr << " leftDown="  << blob.toSymbol(leftDown)  << " leftAccum="  << leftAccum  << " leftAllStates="  << blob.toSymbol(leftAllStates);
//                  std::cerr << " rightDown=" << blob.toSymbol(rightDown) << " rightAccum=" << rightAccum << " rightAllStates=" << blob.toSymbol(rightAllStates) << " \n";
                    addToAncProbVec(*ancVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, instance);
                }
            }
        }
    }
    return probsAdded;
}
*/

/*
void ProbInfo::calculateSymmetric(const ProbInfo & leftPI, double leftEdgeLen, const ProbInfo & rightPI, double rightEdgeLen,
                                  TiMatFunc fn, const MTInstance & instance)
{
  fn(leftEdgeLen, GetPatData(0).firstMatVec.GetAlias());
  const double *** leftPMatVec = const_cast<const double ***>(GetPatData(0).firstMatVec.GetAlias());
  fn(rightEdgeLen, GetPatData(0).secMatVec.GetAlias());
  const double *** rightPMatVec = const_cast<const double ***>(GetPatData(0).secMatVec.GetAlias());

  const unsigned int leftMaxP = leftPI.getMaxParsScore();
  const unsigned int rightMaxP = rightPI.getMaxParsScore();
  const unsigned int maxparscore = 1 + leftMaxP + rightMaxP;

  unsigned int nRat = GetPatData(0).GetNumRates();
  unsigned int nStat = GetPatData(0).GetNumStates();

  this->byParsScore.clear();
  this->byParsScore.resize(maxparscore + 1);

  this->nLeavesBelow = leftPI.getNLeavesBelow() + rightPI.getNLeavesBelow();

  unsigned obsmaxparscore = 0;
  if (true) {
    const ProbForParsScore & leftFPS = leftPI.getByParsScore(0);
        const ProbForParsScore & rightFPS = rightPI.getByParsScore(0);
        ProbForParsScore & forZeroSteps = this->byParsScore[0];

        unsigned stateIndex = 0;
        const BitField sc = 1;

        const std::vector<double> * leftProbs = leftFPS.getProbsForDownPassAndObsMask(sc, sc);
        assert(leftProbs != 0L);
        const std::vector<double> * rightProbs = rightFPS.getProbsForDownPassAndObsMask(sc, sc);
        assert(rightProbs != 0L);

        std::vector<double> & toFillVec = forZeroSteps.byDownPass[sc][sc];
        toFillVec.assign(nRat*nStat, 0.0);

        addToAncProbVec(toFillVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, instance);
    }
  unsigned currscore = 1;
  std::cerr << "from line:" << __LINE__ << ": currscore = " << currscore << ".\n";
  bool scObserved = true;
  ProbForParsScore & forCurrScore = this->byParsScore[currscore];
  BitField downPass = 1;
  if (true) {
    unsigned int leftAccum = 0;
        const ProbForParsScore * leftFPS = &(leftPI.getByParsScore(leftAccum));
        unsigned int rightAccum = 1; // one step on the right side
        const ProbForParsScore * rightFPS = &(rightPI.getByParsScore(rightAccum));

        for (int leftright = 0; leftright < 2; ++leftright) {
            BitField leftDown = 1; // set {0}  is 1 in our bitfield notation
            const MaskToProbsByState * leftM2PBS = leftFPS->getMapPtrForDownPass(leftDown);
            assert(leftM2PBS != 0L);
            const BitField leftAllStates = 1; // set {0}  is 1 in our bitfield notation
            const std::vector<double> * leftProbs = getProbsForStatesMask(leftM2PBS, leftAllStates);
            assert(leftProbs != 0L);
            const BitField rightAllStates = (1|2); // set {0}  is 1 in our bitfield notation
            BitField rightDownArray[] = {1, 1|2};
            for (int rdi = 0 ; rdi < 2; ++rdi) {
                BitField rightDown = rightDownArray[rdi];
                const MaskToProbsByState * rightM2PBS = rightFPS->getMapPtrForDownPass(rightDown);
                assert(rightM2PBS != 0L);
                const std::vector<double> * rightProbs = getProbsForStatesMask(rightM2PBS, rightAllStates);
                assert(rightProbs != 0L);

                const BitField ancAllField = (1|2);
                MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];
                std::vector<double> * ancVec = getMutableProbsForStatesMask(&forCurrScoreDownPass, ancAllField);
                if (ancVec == 0L) { // if we have not visited this set of probabilities, start with a vector of 0's
                    ancVec = &(forCurrScoreDownPass[ancAllField]); // get the memory
                    ancVec->assign(nRat*nStat, 0.0); // set it to 0.0
                }
                addToAncProbVec(*ancVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, instance);
            }
            std::swap(leftFPS, rightFPS);
        }
  }
  // if the ancestor has state set of {0,1}, and currScore = 1 then
    //  both children must display a single state.
    downPass = 3;
    unsigned int leftAccum = 0;
    BitField leftDown = 1; // set {0}  is 1 in our bitfield notation
    unsigned int rightAccum = 0;
    BitField rightDown = 1; // set {0}  is 1 in our bitfield notation
    const ProbForParsScore & leftFPS = leftPI.getByParsScore(leftAccum);
    const MaskToProbsByState * leftM2PBS = leftFPS.getMapPtrForDownPass(leftDown);
    assert(leftM2PBS != 0L);
    const ProbForParsScore & rightFPS = rightPI.getByParsScore(rightAccum);
    const MaskToProbsByState * rightM2PBS = rightFPS.getMapPtrForDownPass(rightDown);
    assert(rightM2PBS != 0L);
    const BitField leftAllStates = 1; // set {0}  is 1 in our bitfield notation
    const BitField rightAllStates = 1; // set {0}  is 1 in our bitfield notation
    const std::vector<double> * leftProbs = getProbsForStatesMask(leftM2PBS, leftAllStates);
    assert(leftProbs != 0L);
    const std::vector<double> * rightProbs = getProbsForStatesMask(rightM2PBS, rightAllStates);
    assert(rightProbs != 0L);

  // even though, we are accessing state set {0} for each child, we are using
    // symmetry to treat one of the children as having state {1}.
    // thus, our ancestor has observed state set {0, 1} (or 3 in BitField notation).
    const BitField ancAllField = 3;
    MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[0];
    std::vector<double> * ancVec = getMutableProbsForStatesMask(&forCurrScoreDownPass, ancAllField);
    if (ancVec == 0L) { // if we have not visited this set of probabilities, start with a vector of 0's
        ancVec = &(forCurrScoreDownPass[ancAllField]); // get the memory
        ancVec->assign(nRat*nStat, 0.0); // set it to 0.0
    }
    std::vector<unsigned int> stateCodeTranslationVec;
    for (unsigned i = 0; i < nStat; ++i) {
        stateCodeTranslationVec.push_back(i);
    }
    stateCodeTranslationVec[0] = 1;
  stateCodeTranslationVec[1] = 0;
    addToAncProbVecSymmetric(*ancVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, stateCodeTranslationVec, instance);
    addToAncProbVecSymmetric(*ancVec, rightPMatVec, rightProbs, leftPMatVec, leftProbs, stateCodeTranslationVec, instance);
*/
  /*
    for (BitField downPass = 1; ; ++downPass) {
        const unsigned numStatesInMask = blob.getNumStates(downPass);
        if (numStatesInMask - 1 <= currScore) { // we cannot demand 3 states seen, but only 1 parsimony change... (all the probs will be zero, so we can skip them)

            MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];

            if (blob.getNumStates(downPass) > 1) {
                std::cerr << "from line: " << __LINE__<< ": downPass = " << (int)downPass << " " << blob.toSymbol(downPass) << " UNIONS:\n";
                const VecMaskPair & forUnions = blob.pairsForUnionForEachDownPass[downPass];
                bool didCalculations = this->allCalcsForAllPairs(forCurrScoreDownPass,
                                          forUnions,
                                          leftPI,
                                          leftPMatVec,
                                          rightPI,
                                          rightPMatVec,
                                          currScore - 1,
                                          false,
                                          blob);
                if (didCalculations)
                    scObserved = true;
            }
            if (leftMaxP + rightMaxP >= currScore) {
                std::cerr << "from line: " << __LINE__<< ": downPass = " << blob.toSymbol(downPass) << " INTERSECTIONS:\n";
                const VecMaskPair & forIntersections = blob.pairsForIntersectionForEachDownPass[downPass];
                bool didCalculations = this->allCalcsForAllPairs(forCurrScoreDownPass,
                                          forIntersections,
                                          leftPI,
                                          leftPMatVec,
                                          rightPI,
                                          rightPMatVec,
                                          currScore,
                                          true,
                                          blob) || scObserved;
                if (didCalculations)
                    scObserved = true;
            }
        }
        if (downPass == blob.lastBitField)
            break;
        assert(downPass < blob.lastBitField);
    }
    */
  /*
    if (scObserved)
        obsmaxparscore = currscore;

    for (currscore = 2;  currscore <= maxparscore; ++currscore) {
        std::cerr << "from line: " << __LINE__<< ": currscore = " << currscore << ".\n";
        bool scObserved = false;
        ProbForParsScore & forCurrScore = this->byParsScore[currscore];
        for (BitField downPass = 1; ; ++downPass) {
            const unsigned numStatesInMask = GetPatData(0).getnstates(downPass);
            if (numStatesInMask - 1 <= currscore) { // we cannot demand 3 states seen, but only 1 parsimony change... (all the probs will be zero, so we can skip them)

                MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];

                if (GetPatData(0).getnstates(downPass) > 1) {
                    std::cerr << "from line: " << __LINE__<< ": downPass = " << (int)downPass << " " << GetPatData(0).toSymbol(downPass) << " UNIONS:\n";
                    const VecMaskPair & forUnions = GetPatData(0).pairsForUnionForEachDownPass[downPass];
                    scObserved = this->allCalcsForAllPairs(forCurrScoreDownPass,
                                              forUnions,
                                              leftPI,
                                              leftPMatVec,
                                              rightPI,
                                              rightPMatVec,
                                              currscore - 1,
                                              false,
                                              instance) || scObserved;
                }
                if (leftMaxP + rightMaxP >= currscore) {
                    std::cerr << "from line: " << __LINE__<< ": downPass = " << GetPatData(0).toSymbol(downPass) << " INTERSECTIONS:\n";
                    const VecMaskPair & forIntersections = GetPatData(0).pairsForIntersectionForEachDownPass[downPass];
                    scObserved = this->allCalcsForAllPairs(forCurrScoreDownPass,
                                              forIntersections,
                                              leftPI,
                                              leftPMatVec,
                                              rightPI,
                                              rightPMatVec,
                                              currscore,
                                              true,
                                              instance) || scObserved;
                }
            }
            if (downPass == GetPatData(0).lastBitField)
                break;
            assert(downPass < GetPatData(0).lastBitField);
        }
        if (scObserved)
            obsmaxparscore = currscore;
    }
    if (obsmaxparscore < maxparscore) {
        std::cerr << "from line: " << __LINE__<< '\n'; std::cerr << "maxParsScore  = " << maxparscore << " obsMaxParsScore = " << obsmaxparscore << "\n";
        this->byParsScore.resize(obsmaxparscore + 1);
    }
}
*/

/*
void ProbInfo::calculate(const ProbInfo & leftPI, double leftEdgeLen,
                                   const ProbInfo & rightPI, double rightEdgeLen,
                                   TiMatFunc fn, const MTInstance & instance) {

  fn(leftEdgeLen, GetPatData(0).firstMatVec.GetAlias());
  const double *** leftPMatVec = const_cast<const double ***>(GetPatData(0).firstMatVec.GetAlias());
  fn(rightEdgeLen, GetPatData(0).secMatVec.GetAlias());
  const double *** rightPMatVec = const_cast<const double ***>(GetPatData(0).secMatVec.GetAlias());

  const unsigned int leftMaxP = leftPI.getMaxParsScore();
  const unsigned rightMaxP = rightPI.getMaxParsScore();
  const unsigned maxParsScore = 1 + leftMaxP + rightMaxP;
  unsigned int nRat = GetPatData(0).GetNumRates();
  unsigned int nStat = GetPatData(0).GetNumStates();

  this->byParsScore.clear();
  this->byParsScore.resize(maxParsScore + 1); // add one to account zero
  this->nLeavesBelow = leftPI.getNLeavesBelow() + rightPI.getNLeavesBelow();
  unsigned obsMaxParsScore = 0;
  // do the calculations for staying in the constant patterns, these are more simple than the general calcs...
  if (true) { //@TEMP if true so that variables are scoped.
    const ProbForParsScore & leftFPS = leftPI.getByParsScore(0);
    const ProbForParsScore & rightFPS = rightPI.getByParsScore(0);
    ProbForParsScore & forZeroSteps = this->byParsScore[0];

    unsigned stateIndex = 0;
    for (std::vector<BitField>::const_iterator scIt = GetPatData(0).singleStateCodes.begin();
         scIt != GetPatData(0).singleStateCodes.end();
         ++scIt, ++stateIndex) {
      const BitField sc = *scIt;
      const std::vector<double> * leftProbs = leftFPS.getProbsForDownPassAndObsMask(sc, sc);
      assert(leftProbs != 0L);
      const std::vector<double> * rightProbs = rightFPS.getProbsForDownPassAndObsMask(sc, sc);
      assert(rightProbs != 0L);

      std::vector<double> & toFillVec = forZeroSteps.byDownPass[sc][sc];
      toFillVec.assign(nRat*nStat, 0.0);

      addToAncProbVec(toFillVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, instance);
    }
  }

  // order N
  for (unsigned currScore = 1; currScore <= maxParsScore; ++currScore) {
    std::cerr << "from line: " << __LINE__<< ": currScore = " << currScore << ".\n";
    bool scObserved = false;
    ProbForParsScore & forCurrScore = this->byParsScore[currScore];
    for (BitField downPass = 1; ; ++downPass) {
      const unsigned numStatesInMask = GetPatData(0).getnstates(downPass);
      if (numStatesInMask - 1 <= currScore) { // we cannot demand 3 states seen, but only 1 parsimony change... (all the probs will be zero, so we can skip them)

        MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];

        if (GetPatData(0).getnstates(downPass) > 1) {
          std::cerr << "from line: " << __LINE__<< ": downPass = " << (int)downPass << " " << GetPatData(0).toSymbol(downPass) << " UNIONS:\n";
          const VecMaskPair & forUnions = GetPatData(0).pairsForUnionForEachDownPass[downPass];
          scObserved = this->allCalcsForAllPairs(forCurrScoreDownPass,
                                                                      forUnions,
                                                                      leftPI,
                                                                      leftPMatVec,
                                                                      rightPI,
                                                                      rightPMatVec,
                                                                      currScore - 1,
                                                                      false,
                                                                      instance) || scObserved;
        }
        if (leftMaxP + rightMaxP >= currScore) {
          std::cerr << "from line: " << __LINE__<< ": downPass = " << GetPatData(0).toSymbol(downPass) << " INTERSECTIONS:\n";
          const VecMaskPair & forIntersections = GetPatData(0).pairsForIntersectionForEachDownPass[downPass];
          scObserved = this->allCalcsForAllPairs(forCurrScoreDownPass,
                                                                      forIntersections,
                                                                      leftPI,
                                                                      leftPMatVec,
                                                                      rightPI,
                                                                      rightPMatVec,
                                                                      currScore,
                                                                      true,
                                                                      instance) || scObserved;
        }
      }
      if (downPass == GetPatData(0).lastBitField)
        break;
      assert(downPass < GetPatData(0).lastBitField);
    }
    if (scObserved)
      obsMaxParsScore = currScore;
  }
  if (obsMaxParsScore < maxParsScore) {
    std::cerr << "from line: " << __LINE__<< '\n'; std::cerr << "maxParsScore  = " << maxParsScore << " obsMaxParsScore = " << obsMaxParsScore << "\n";
    this->byParsScore.resize(obsMaxParsScore + 1);
  }

}
*/

/*
unsigned PatternSummary::incrementCount(unsigned s, BitField m, unsigned toAdd) {
    if (s >= this->byParsScore.size())
        this->byParsScore.resize(s + 1);
    BitsToCount & mapper = this->byParsScore[s];
    BitsToCount::iterator mIt = mapper.find(m);
    if (mIt == mapper.end()) {
        mapper[m] = toAdd;
        return toAdd;
    }
    unsigned prev = mIt->second;
    mIt->second = toAdd + prev;
    return toAdd + prev;
}
*/
/*
ExpectedPatternSummary::ExpectedPatternSummary(const ProbInfo & rootProbInfo, const MTInstance & instance) {
    if (instance.GetCharModel().isMkvSymm) {
        const unsigned maxNumSteps = rootProbInfo.getMaxParsScore();
        std::vector<double> emptyRow(instance.GetCharModel().lastBitField + 1, 0.0);
        this->probsByStepsThenObsStates.resize(maxNumSteps + 1, emptyRow);

        //For constant patterns
        const ProbForParsScore & constFPS = rootProbInfo.getByParsScore(0);
        const std::vector<double> * pVec = constFPS.getProbsForDownPassAndObsMask(1,1);
        assert(pVec);
        double patClassProb = 0.0;
        std::vector<double>::const_iterator wtIt = GetPatData(0).categStateProb.begin();
        std::vector<double>::const_iterator pIt = pVec->begin();
        for (; wtIt != GetPatData(0).categStateProb.end(); ++wtIt, ++pIt) {
            assert(pIt != pVec->end());
            patClassProb += (*wtIt) * (*pIt);
        }
        for (std::vector<BitField>::const_iterator scIt = GetPatData(0).singleStateCodes.begin();
            scIt != GetPatData(0).singleStateCodes.end(); ++scIt) {
            const BitField downPass = *scIt;
            this->probsByStepsThenObsStates[0][downPass] = patClassProb;
            }

        //Other Patterns
        for (unsigned i = 1; i <= maxNumSteps; ++i) {
            const ProbForParsScore & fps = rootProbInfo.getByParsScore(i);
            for (BitField downPass = 1; ; ++downPass) {
                for (BitField obsStates = 1; ; ++obsStates) {
                    double patClassProb = 0.0;
                    std::vector<double>::const_iterator wtIt = GetPatData(0).categStateProb.begin();
                    const std::vector<double> * pVec = fps.getProbsForDownPassAndObsMask(downPass, obsStates);
                    if (pVec) {
                        std::vector<double>::const_iterator pIt = pVec->begin();
                        for (; wtIt != GetPatData(0).categStateProb.end(); ++wtIt, ++pIt) {
                            assert(pIt != pVec->end());
                            patClassProb += (*wtIt) * (*pIt);
                        }
                    }
                    this->probsByStepsThenObsStates[i][obsStates] += patClassProb;
                    if (obsStates == GetPatData(0).lastBitField)
                        break;
                }
                if (downPass == GetPatData(0).lastBitField)
                    break;
            }
        }
    } else {
        const unsigned maxNumSteps = rootProbInfo.getMaxParsScore();
        std::vector<double> emptyRow(GetPatData(0).lastBitField + 1.0, 0.0);
        this->probsByStepsThenObsStates.resize(maxNumSteps + 1, emptyRow);
        const ProbForParsScore & constFPS = rootProbInfo.getByParsScore(0);
        for (std::vector<BitField>::const_iterator scIt = GetPatData(0).singleStateCodes.begin();
            scIt != GetPatData(0).singleStateCodes.end();
            ++scIt) {
            const BitField downPass = *scIt;
            const std::vector<double> * pVec = constFPS.getProbsForDownPassAndObsMask(downPass, downPass);
            assert(pVec);
            double patClassProb = 0.0;
            std::vector<double>::const_iterator wtIt = GetPatData(0).categStateProb.begin();
            std::vector<double>::const_iterator pIt = pVec->begin();
            for (; wtIt != GetPatData(0).categStateProb.end(); ++wtIt, ++pIt) {
                assert(pIt != pVec->end());
                patClassProb += (*wtIt) * (*pIt);
            }
            this->probsByStepsThenObsStates[0][downPass] = patClassProb;
            }
            for (unsigned i = 1; i <= maxNumSteps; ++i) {
                const ProbForParsScore & fps = rootProbInfo.getByParsScore(i);
                for (BitField downPass = 1; ; ++downPass) {
                    for (BitField obsStates = 1; ; ++obsStates) {
                        double patClassProb = 0.0;
                        std::vector<double>::const_iterator wtIt = GetPatData(0).categStateProb.begin();
                        const std::vector<double> * pVec = fps.getProbsForDownPassAndObsMask(downPass, obsStates);
                        if (pVec) {
                            std::vector<double>::const_iterator pIt = pVec->begin();
                            for (; wtIt != GetPatData(0).categStateProb.end(); ++wtIt, ++pIt) {
                                assert(pIt != pVec->end());
                                patClassProb += (*wtIt) * (*pIt);
                            }
                        }
                        this->probsByStepsThenObsStates[i][obsStates] += patClassProb;
                        if (obsStates == GetPatData(0).lastBitField)
                            break;
                    }
                    if (downPass == GetPatData(0).lastBitField)
                        break;
                }
            }
    }
}
*/

// Maybe to replace ProbInfo or other data structure
// Class to store node data for uninformative case
// NodeDataStructure in Uninformative_case.cpp
class NodeInfo {
    public:
        NodeInfo(unsigned numStates) {
          numLeaves = 0;
          int len = 1 << numStates; // bitwise shift because character states as bitfields
          ProbForObsStateSet dummy (numStates);
          probVec.assign(len, dummy);
        }

        ProbForObsStateSet & getForObsStateSet(int obs) {
            return probVec.at(obs);
        }

        int getNumLeaves() {
          return numLeaves;
        }

        void setNumLeaves(int n) {
          numLeaves = n;
        }
        
        void setIsMissing(bool val) {
          isMissing = val;
        }
        bool getIsMissing() const {
            return isMissing;
        }
    private:
        std::vector<ProbForObsStateSet> probVec;
        int numLeaves;
        bool isMissing;
};

int countBits(int x)
{
    int num = 0;
    while(x > 0)
    {
        if(x& 1)
            ++num;
        x = (x>>1);
    }
    return num;
}

int convertIndexToBit(int ind) {
  return 1 << ind;
}

int convertBitToIndex(int i) {
  int ind = 0;
  while(i > 0) {
    if(i == 1)
      return ind;
    if(1& i) {
      std::cerr << "Illegal bit value \n";
      exit(1);
    }
    ind++;
    i = (i>>1);
  }
  std::cerr << "Zero to Convert Bits \n";
  exit(1);
}

std::vector<int> subsetsContainingGivenState(int fullSet, int givenState) {
  std::set<int> subsets;
  int i = 1;
  while(i <= fullSet) {
    int j = i& fullSet;
    if(j& givenState)
      subsets.insert(j);
    i++;
  }
  return std::vector<int> (subsets.begin(), subsets.end());
}

std::vector<int> subsetsOfGivenSize(int obsStSet, int numBits) {
  std::set<int> subsets;
  int i = 1;
  while(i <= obsStSet) {
    int j = 1& obsStSet;
    if(countBits(j) == numBits)
      subsets.insert(j);
    i++;
  }
  return std::vector<int> (subsets.begin(), subsets.end());
}

int getNextCommStSet(const int obsStSet, int i) {
  int ind, binRep;
  if(i == 1) {
    ind = 0;
    binRep = 1;
  } else {
    ind = i + 1;
    binRep = 1 << ind;
    if(binRep > obsStSet)
      return -2;
  }
  while((binRep & obsStSet) == 0) {
    binRep <<= 1;
    ind++;
  }
  return ind;
}

double pclassCalcTransitionProb(int ancIndex, int i, double edgeLen, MTInstance & instance){
  double * tiVec = GetPatData(0).calcTransitionProb(edgeLen);
  int nStates = GetPatData(0).GetNumStates();
  return tiVec[ancIndex*nStates + i];
}

double calcProbOfSubtreeForObsStSetAndComm(NodeInfo * subtreeInfo, int ancIndex, int obsBits, int commonStates, double edgeLen, MTInstance &instance) {
  double p = 0.0;
  ProbForObsStateSet & childProbSet = subtreeInfo->getForObsStateSet(obsBits);
  std::vector<double> & childProb = childProbSet.getProbForCommState(commonStates);
  for(int i = 0; i < GetPatData(0).GetNumStates(); i++)  {
    double transProb = pclassCalcTransitionProb(ancIndex, i, edgeLen, instance);
    double partialLike = childProb[i];
    double x = transProb * partialLike;
    p += x;
  }
  return p;
}

double calcProbOfSubtreeForObsStSetNoRepeated(NodeInfo * subtreeInfo, int ancIndex, int obsBits, double edgeLen, MTInstance &instance){
  return calcProbOfSubtreeForObsStSetAndComm(subtreeInfo, ancIndex, obsBits, -1, edgeLen, instance);
}

// Traverse the tree (postorder) and calculate pattern class probabilities
// might not need this
/*
void calcPatternClassProbs(MTInstance &instance, TiMatFunc fn)
{
    initInfo(instance);
    ProbInfo * rootpinfo;
    bool needToDelRootProbInfo = false;
    ProbInfo tiprobinfo;
    NodeIDToProbInfo nodeIDToProbInfo;

    Node * vRoot = instance.tree.GetRoot();
    vRoot = vRoot->leftChild->rightSib;
    PostorderForNodeIterator postTravIter = postorder(vRoot);
    Arc postTravArc = postTravIter.get();
    tiprobinfo.createForTip(instance);
    try {
      while (postTravArc.toNode) {
        const Node * nd = postTravArc.fromNode;
        std::vector<Node *> children = nd->GetChildren();
        const unsigned numChildren = children.size();
        NodeID currNdID(nd, 0);
        if (numChildren == 0){
          nodeIDToProbInfo[currNdID] = &tiprobinfo;
        }
        else {
          if (numChildren == 1){
            std::cout << "Trees of degree 2 are not supported\n";
            throw;
          }
          ProbInfo * currProbInfo = new ProbInfo();
          nodeIDToProbInfo[currNdID] = currProbInfo;
          const Node * leftNd = children[0];
          const Node * rightNd = children[1];
          NodeIDToProbInfo::const_iterator leftPIIt= nodeIDToProbInfo.find(NodeID(leftNd, 0));
          assert(leftPIIt != nodeIDToProbInfo.end());
          ProbInfo * leftPI = leftPIIt->second;
          assert(leftPI);
          NodeIDToProbInfo::const_iterator rightPIIt= nodeIDToProbInfo.find(NodeID(rightNd, 0));
          assert(rightPIIt != nodeIDToProbInfo.end());
          ProbInfo * rightPI = rightPIIt->second;
          assert(rightPI != 0L);
          ProbInfo lt, rt;
          if (GetPatData(0).isMkvSymm) {
            currProbInfo->calculateSymmetric(*leftPI, leftNd->GetEdgeLen(), *rightPI, rightNd->GetEdgeLen(), fn, instance);
          }
          else {
            currProbInfo->calculate(*leftPI, leftNd->GetEdgeLen(), *rightPI, rightNd->GetEdgeLen(), fn, instance);
          }
          if (leftPI->getNLeavesBelow() > 1) {
            delete leftPI;
            nodeIDToProbInfo[NodeID(leftNd, 0)] = 0L;
          }
          if (rightPI->getNLeavesBelow() > 1) {
            delete rightPI;
            nodeIDToProbInfo[NodeID(rightNd, 0)] = 0L;
          }
          if (numChildren > 2) {
            if (nd != instance.tree.GetRoot() || numChildren > 3) {
              std::cout << "Parsimony scoring on non-binary trees is not supported\n";
              throw;
            }
            const Node * lastNd = children.at(2);
            NodeIDToProbInfo::const_iterator lastPIIt = nodeIDToProbInfo.find(NodeID(lastNd, 0));
            assert(lastPIIt != nodeIDToProbInfo.end());
            ProbInfo * lastPI = lastPIIt->second;
            assert(lastPI);
            rootpinfo = new ProbInfo();
            needToDelRootProbInfo = true;
            if (GetPatData(0).isMkvSymm) {
              rootpinfo->calculateSymmetric(*currProbInfo, 0.0, *lastPI, lastNd->GetEdgeLen(), fn, instance);
            }
            else {
              rootpinfo->calculate(*currProbInfo, 0.0, *lastPI, lastNd->GetEdgeLen(), fn, instance);
            }
          } else if (nd == instance.tree.GetRoot())
            rootpinfo = currProbInfo;
        }
        _DEBUG_VAL(postTravArc.fromNode->number);
        // advance arc
        postTravArc = postTravIter.next();
      }
      assert(rootpinfo != 0L);
      //const ExpectedPatternSummary eps(*rootpinfo, instance);
      //eps.write(std::cout, instance);
    }
    catch (...) {
      PostorderForNodeIterator pInfoFreer1 = postorder(vRoot);
      freeProbInfo(pInfoFreer1, nodeIDToProbInfo);
      if (needToDelRootProbInfo)
        delete rootpinfo;
      throw;
    }
    PostorderForNodeIterator pInfoFreer2 = postorder(vRoot);
    freeProbInfo(pInfoFreer2, nodeIDToProbInfo);
    if (needToDelRootProbInfo)
      delete rootpinfo;
}
*/

// calculate probabilities of uninformative patterns
void calcUninformativePatterns(MTInstance & instance)
{
  Node * nd = instance.tree.GetRoot();
  PostorderForNodeIterator poTrav = postorder(nd);
  Arc arc = poTrav.get();
  std::map<Node *, NodeInfo *> nodeToInfoMap;
  unsigned numStates = GetPatData(0).GetNumStates();
  NodeInfo * currNdInfo = 0L;
  assert(arc.toNode);
  while(arc.toNode) {
    Node * currNd = arc.fromNode;
    std::vector<Node *> children = currNd->GetChildren();
    const unsigned numChildren = children.size();
    currNdInfo = new NodeInfo(numStates);
    NodeID currNdID(currNd, 0);
    nodeToInfoMap[currNd] = currNdInfo;
    if (numChildren == 0) {
      for(int i = 0; i < numStates; i++) {
        int ss=1 << i;
        ProbForObsStateSet & p = currNdInfo->getForObsStateSet(ss);
        std::vector<double> & v = p.getProbForCommState(-1);
        v[i] = 1.0;
        currNdInfo->setNumLeaves(1);
      }
    } else {

    if (numChildren != 2) {
      std::cerr << "Trees must be binary \n";
      exit(1);
    }

      Node * leftChild = children[0];
      NodeInfo * leftNdInfo = nodeToInfoMap[leftChild];
      Node * rightChild = children[1];
      NodeInfo * rightNdInfo = nodeToInfoMap[rightChild];
      if(leftNdInfo->getIsMissing() || rightNdInfo->getIsMissing()) {
        if(leftNdInfo->getIsMissing() && rightNdInfo->getIsMissing()) {
          currNdInfo->setIsMissing(true);
        } else {
          // add edge lengths and progress to next node
        }
      } else {}
      currNdInfo->setNumLeaves(leftNdInfo->getNumLeaves() + rightNdInfo->getNumLeaves());

      stateSetContainer::const_iterator ssCit = GetPatData(0).stateSetBegin();
      for (; ssCit != GetPatData(0).stateSetEnd(); ssCit++) {
        const int & obsStSet = *ssCit;
        int common = -1;   // this is 11... in bits
        int numObsSt = countBits(obsStSet);

        while(common>-2) {
          ProbForObsStateSet & currNdProbSet = currNdInfo->getForObsStateSet(obsStSet);
          std::vector<double> & currNdProbVec = currNdProbSet.getProbForCommState(common);

          if(common == -1) {
            if (currNdInfo->getNumLeaves() == numObsSt) {
              for (int anc = 0; anc < numStates; anc++) {
                std::cerr << "ObsStSet " << obsStSet << '\n';
                currNdProbVec[anc] = 0.0;
                std::vector<int> leftObsStSets = subsetsOfGivenSize(obsStSet, leftNdInfo->getNumLeaves());
                for (int j = 0; j < leftObsStSets.size(); j++) {
                  int leftObsStSet = leftObsStSets[j];
                  int rightObsStSet = obsStSet - leftObsStSet;

                  double leftProb, rightProb;
                  double leftedgeLen = leftChild->GetEdgeLen();
                  if(leftNdInfo->getNumLeaves() == 1) {
                    leftProb = pclassCalcTransitionProb(anc, convertBitToIndex(leftObsStSet), leftedgeLen, instance);
                  } else {
                    leftProb = calcProbOfSubtreeForObsStSetNoRepeated(leftNdInfo, anc, leftObsStSet, leftedgeLen, instance);
                  }

                  double rightEdgeLen = rightChild->GetEdgeLen();
                  if(rightNdInfo->getNumLeaves() == 1) {
                    rightProb = pclassCalcTransitionProb(anc, convertBitToIndex(rightObsStSet), rightEdgeLen, instance);
                  } else {
                    rightProb = calcProbOfSubtreeForObsStSetNoRepeated(rightNdInfo, anc, rightObsStSet, rightEdgeLen, instance);
                  }
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;
                  }
                }
              }
            } else {

              int commonBits = convertIndexToBit(common);
              for (int anc = 0; anc < numStates; anc++) {

                currNdProbVec[anc] = 0.0;
                int leftCommSt, rightCommSt;
                leftCommSt = common;
                rightCommSt = common;
                std::vector<int> obsStSetsWithComm = subsetsContainingGivenState(obsStSet, commonBits);

                for(int j=0; j < obsStSetsWithComm.size(); j++) {
                  int leftObsStSet = obsStSetsWithComm[j];
                  int rightObsStSet = obsStSet - leftObsStSet + commonBits;

                  double leftProb, rightProb;
                  double leftEdgeLen = leftChild->GetEdgeLen();
                  leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, common, leftEdgeLen, instance);
                  double rightEdgeLen = rightChild->GetEdgeLen();
                  rightEdgeLen = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, common, rightEdgeLen, instance);
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;
                }

                leftCommSt = -1;
                rightCommSt = common;
                //add probability when only right common, left not repeated
                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                  int rightObsStSet = obsStSetsWithComm[j];
                  int leftObsStSet = obsStSet - rightObsStSet + commonBits;
                  double leftProb, rightProb;
                  double leftEdgeLen = leftChild->GetEdgeLen();
                  leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, -1, leftEdgeLen, instance);
                  double rightEdgeLen = rightChild->GetEdgeLen();
                  rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, common, rightEdgeLen, instance);
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;

                  //Now consider when the left is not displayed by commonBits as observed States
                  leftObsStSet = obsStSet - rightObsStSet;
                  if(leftObsStSet != 0) {
                    leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, -1, leftEdgeLen, instance);
                    jointNdProb = leftProb * rightProb;
                    currNdProbVec[anc] += jointNdProb;
                    }
                }

                leftCommSt = common;
                rightCommSt = -1;
                //add probability when only left common
                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                  int rightObsStSet = obsStSetsWithComm[j];
                  int leftObsStSet = obsStSet - rightObsStSet + commonBits;
                  double leftProb, rightProb;
                  double leftEdgeLen = leftChild->GetEdgeLen();
                  leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, common, leftEdgeLen, instance);
                  double rightEdgeLen = rightChild->GetEdgeLen();
                  rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, -1, rightEdgeLen, instance);
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;

                  //Now consider when the right is not displayed by commonBits as observed States
                  rightObsStSet = obsStSet - leftObsStSet;
                  if(rightObsStSet != 0) {
                    rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, -1, rightEdgeLen, instance);
                    jointNdProb = leftProb * rightProb;
                    currNdProbVec[anc] += jointNdProb;
                  }
                }

                leftCommSt = -1;
                rightCommSt = -1;
                //add probability when neither common
                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                  int rightObsStSet = obsStSetsWithComm[j];
                  int leftObsStSet = obsStSet - rightObsStSet + commonBits;
                  double leftProb, rightProb;
                  double leftEdgeLen = leftChild->GetEdgeLen();
                  leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, -1, leftEdgeLen, instance);
                  double rightEdgeLen = rightChild->GetEdgeLen();
                  rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, -1, rightEdgeLen, instance);
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;
                }
              }
            }
            common = getNextCommStSet(obsStSet, common);
          }
        }
      }
      arc = poTrav.next();
    }
    // return currNdInfo
    //do something else

  } // calcUninformativePatterns function end


} // namespace
