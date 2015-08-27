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

namespace mt {

//MACROS
#define GetPatData      instance.GetCharModel()

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
	for (std::vector<BitField>::const_iterator scIt = GetPatData.singleStateCodes.begin();
		 	scIt != GetPatData.singleStateCodes.end();
		 	++scIt, ++stateIndex) {
		const BitField sc = *scIt;
		std::vector<double> & pVec = forZeroSteps.byDownPass[sc][sc];
		pVec.assign(GetPatData.GetNumRates()*GetPatData.GetNumStates(), 0.0);
		for (unsigned r = 0; r < GetPatData.GetNumRates(); ++r)
			pVec[GetPatData.GetNumStates()*r + stateIndex] = 1.0;
	}
	this->nLeavesBelow = 1;
}

void ProbInfo::addToAncProbVecSymmetric(std::vector<double> & pVec,
		const double *** leftPMatVec, const std::vector<double> * leftProbs,
		const double *** rightPMatVec, const std::vector<double> * rightProbs,
		const std::vector<unsigned int> & rightChildStateCodeTranslation,
		const MTInstance & instance) {
	if (leftProbs == 0L || rightProbs == 0L)
		return;
	unsigned rOffset = 0;
	// ignore loop over rates blob.nRates = 1
	for (unsigned r = 0; r < GetPatData.GetNumRates(); ++r) {
		// this code looks up the correct transition prob matrix
		const double ** leftPMat = leftPMatVec[r];
		const double ** rightPMat = rightPMatVec[r];


		for (unsigned ancState = 0; ancState < GetPatData.GetNumStates(); ++ancState) {
			double leftProb = 0.0;
			double rightProb = 0.0;
			for (unsigned desState = 0; desState < GetPatData.GetNumStates(); ++desState) {

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
		rOffset += GetPatData.GetNumStates();
	}
}

void ProbInfo::addToAncProbVec(std::vector<double> & pVec,
		const double *** leftPMatVec, const std::vector<double> * leftProbs,
		const double *** rightPMatVec, const std::vector<double> * rightProbs,
		const MTInstance & instance) {
	if (leftProbs == 0L || rightProbs == 0L)
		return;
	unsigned rOffset = 0;
	// ignore loop over rates blob.nRates = 1
	for (unsigned r = 0; r < GetPatData.GetNumRates(); ++r) {
		// this code looks up the correct transition prob matrix
		const double ** leftPMat = leftPMatVec[r];
		const double ** rightPMat = rightPMatVec[r];


		for (unsigned ancState = 0; ancState < GetPatData.GetNumStates(); ++ancState) {
			double leftProb = 0.0;
			double rightProb = 0.0;
			for (unsigned desState = 0; desState < GetPatData.GetNumStates(); ++desState) {
				const double leftTiProb = leftPMat[ancState][desState];
				const double leftAccumProb = (*leftProbs)[rOffset + desState];
				leftProb += leftTiProb*leftAccumProb;



				const double rightTiProb = rightPMat[ancState][desState];
				const double rightAccumProb = (*rightProbs)[rOffset + desState];
				rightProb += rightTiProb*rightAccumProb;
			}
			pVec[rOffset + ancState] += leftProb*rightProb;
		}
		rOffset += GetPatData.GetNumStates();
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

void initInfo(MTInstance &instance) {
	GetPatData.isMkvSymm = false;
	GetPatData.pVecLen = GetPatData.GetNumStates()*GetPatData.GetNumRates();
	GetPatData.categStateProb.assign(GetPatData.pVecLen, 1.0/((double)GetPatData.pVecLen));
	GetPatData.singleStateCodes.clear();
	GetPatData.multiStateCodes.clear();
	GetPatData.stateCodesToSymbols.clear();
	unsigned lbfU = (1 << GetPatData.GetNumStates()) - 1;
	GetPatData.stateCodeToNumStates.assign(lbfU + 1, 0);
	GetPatData.lastBitField = BitField(lbfU);
	BitField sc = 1;
	for (;; sc++) {
		const std::set<BitField> sbf = toElements(sc);
		if (sbf.size() == 1) {
			unsigned stInd = GetPatData.singleStateCodes.size();
			GetPatData.singleStateCodes.push_back(sc);
			GetPatData.stateIndexToStateCode.push_back(sc);
			assert(GetPatData.stateIndexToStateCode[stInd] == sc);
			GetPatData.stateCodesToSymbols[sc] = GetPatData.alphabet[stInd]; // FIX THIS
			GetPatData.stateCodeToNumStates.at(sc) = 1;
		} else {
			GetPatData.multiStateCodes.push_back(sc);
			std::string sym;
			for (std::set<BitField>::const_iterator sbfIt = sbf.begin(); sbfIt != sbf.end(); ++sbfIt)
				sym.append(GetPatData.stateCodesToSymbols[*sbfIt]);
			GetPatData.stateCodeToNumStates.at(sc) = sbf.size();
			GetPatData.stateCodesToSymbols[sc] = sym;
		}

		if (sc == GetPatData.lastBitField)
			break;

	}
	GetPatData.pairsForUnionForEachDownPass.clear();
	GetPatData.pairsForUnionForEachDownPass.resize(GetPatData.lastBitField + 1);
	GetPatData.pairsForIntersectionForEachDownPass.clear();
	GetPatData.pairsForIntersectionForEachDownPass.resize(GetPatData.lastBitField + 1);

	for (sc = 1;; sc++) {
		if (GetPatData.stateCodeToNumStates[sc] > 1) {
			VecMaskPair & forUnions = GetPatData.pairsForUnionForEachDownPass[sc];
			for (BitField leftSC = 1; leftSC < sc; ++leftSC) {
				if ((leftSC | sc) != sc)
					continue;
				BitField rightSC = sc - leftSC;
				assert((rightSC | leftSC) == sc);
				forUnions.push_back(MaskPair(leftSC, rightSC));
			}
		}
		VecMaskPair & forIntersections = GetPatData.pairsForIntersectionForEachDownPass[sc];
		for (BitField leftSC = 1; leftSC <= GetPatData.lastBitField ; ++leftSC) {
			for (BitField rightSC = 1; rightSC <= GetPatData.lastBitField ; ++rightSC) {
				if ((leftSC & rightSC) != sc)
					continue;
				forIntersections.push_back(MaskPair(leftSC, rightSC));
			}
		}

		if (sc == GetPatData.lastBitField)
			break;

	}
        GetPatData.statesSupersets.clear();
	GetPatData.statesSupersets.resize(GetPatData.lastBitField + 1);
	for (sc = 1;;++sc) {
		BitFieldRow & ssRow = GetPatData.statesSupersets[sc];
		for (BitField ss = sc;; ++ss) {
			if ((ss & sc) == sc)
				ssRow.push_back(ss);
			if (ss == GetPatData.lastBitField)
				break;
		}

		if (sc == GetPatData.lastBitField)
			break;

	}

	for (sc = 1;; sc++) {
		const VecMaskPair & forUnions = GetPatData.pairsForUnionForEachDownPass[sc];
		const VecMaskPair & forIntersections = GetPatData.pairsForIntersectionForEachDownPass[sc];
		if (sc == GetPatData.lastBitField)
			break;
	}
}


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
{	// order (2^k)^2
	const unsigned leftMaxP = leftPI.getMaxParsScore();
	const unsigned rightMaxP = rightPI.getMaxParsScore();
	bool probsAdded = false;
	for (VecMaskPair::const_iterator fuIt = pairVec.begin(); fuIt != pairVec.end(); ++fuIt) {
		const BitField leftDown = fuIt->first;
		const BitField rightDown = fuIt->second;
		std::cerr << "from line: " << __LINE__<< ": "; std::cerr << "leftDown = " << GetPatData.toSymbol(leftDown) << " rightDown = " << GetPatData.toSymbol(rightDown) << '\n';
		assert(leftDown > 0);
		assert(rightDown > 0);
		if (doingIntersection) {
			assert((leftDown & rightDown) != 0);
		}
		else {
			assert((leftDown & rightDown) == 0);
		}
		const unsigned leftMinAccum = GetPatData.getnstates(leftDown) - 1;
		const unsigned rightMinAccum = GetPatData.getnstates(rightDown) - 1;
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
			const BitFieldRow & leftSSRow = GetPatData.statesSupersets[leftDown];
			// order (2^k)
			for (BitFieldRow::const_iterator lasIt = leftSSRow.begin(); lasIt != leftSSRow.end(); ++lasIt) {
				const BitField leftAllStates = *lasIt;
					std::cerr << "from line: " << __LINE__<< ":	 leftAllStates="  << GetPatData.toSymbol(leftAllStates) << '\n';
				const std::vector<double> * leftProbs = getProbsForStatesMask(leftM2PBS, leftAllStates);
				if (leftProbs == 0L) {
					std::cerr << "from line: " << __LINE__<< ": left child empty bin. Skipping...\n";
					continue;
				}
				const BitFieldRow & rightSSRow = GetPatData.statesSupersets[rightDown];
				// order (2^k)
				for (BitFieldRow::const_iterator rasIt = rightSSRow.begin(); rasIt != rightSSRow.end(); ++rasIt) {
					const BitField rightAllStates = *rasIt;
//					std::cerr << "from line: " << __LINE__<< ":  rightAllStates="  << blob.toSymbol(rightAllStates) << '\n';
					const std::vector<double> * rightProbs = getProbsForStatesMask(rightM2PBS, rightAllStates);
					if (rightProbs == 0L) {
//						std::cerr << "from line: " << __LINE__<< ": right child empty bin. Skipping...\n";
						continue;
					}
					const BitField ancAllField = (rightAllStates|leftAllStates);

					std::vector<double> * ancVec = getMutableProbsForStatesMask(&forCurrScoreDownPass, ancAllField);
					if (ancVec == 0L) {
						ancVec = &(forCurrScoreDownPass[ancAllField]);
						ancVec->assign(GetPatData.GetNumRates()*GetPatData.GetNumStates(), 0.0);
					}
					probsAdded = true;
//					std::cerr << __LINE__ << " adding:";
//					std::cerr << " leftDown="  << blob.toSymbol(leftDown)  << " leftAccum="  << leftAccum  << " leftAllStates="  << blob.toSymbol(leftAllStates);
//					std::cerr << " rightDown=" << blob.toSymbol(rightDown) << " rightAccum=" << rightAccum << " rightAllStates=" << blob.toSymbol(rightAllStates) << " \n";
					addToAncProbVec(*ancVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, instance);
				}
			}
		}
	}
	return probsAdded;
}


void ProbInfo::calculateSymmetric(const ProbInfo & leftPI, double leftEdgeLen, const ProbInfo & rightPI, double rightEdgeLen,
                                  TiMatFunc fn, const MTInstance & instance)
{
  fn(leftEdgeLen, GetPatData.firstMatVec.GetAlias());
  const double *** leftPMatVec = const_cast<const double ***>(GetPatData.firstMatVec.GetAlias());
  fn(rightEdgeLen, GetPatData.secMatVec.GetAlias());
  const double *** rightPMatVec = const_cast<const double ***>(GetPatData.secMatVec.GetAlias());

  const unsigned int leftMaxP = leftPI.getMaxParsScore();
  const unsigned int rightMaxP = rightPI.getMaxParsScore();
  const unsigned int maxparscore = 1 + leftMaxP + rightMaxP;

  unsigned int nRat = GetPatData.GetNumRates();
  unsigned int nStat = GetPatData.GetNumStates();

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
	// 	both children must display a single state.
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
	if (scObserved)
		obsmaxparscore = currscore;

	for (currscore = 2;	 currscore <= maxparscore; ++currscore) {
		std::cerr << "from line: " << __LINE__<< ": currscore = " << currscore << ".\n";
		bool scObserved = false;
		ProbForParsScore & forCurrScore = this->byParsScore[currscore];
		for (BitField downPass = 1; ; ++downPass) {
			const unsigned numStatesInMask = GetPatData.getnstates(downPass);
			if (numStatesInMask - 1 <= currscore) { // we cannot demand 3 states seen, but only 1 parsimony change... (all the probs will be zero, so we can skip them)

				MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];

				if (GetPatData.getnstates(downPass) > 1) {
					std::cerr << "from line: " << __LINE__<< ": downPass = " << (int)downPass << " " << GetPatData.toSymbol(downPass) << " UNIONS:\n";
					const VecMaskPair & forUnions = GetPatData.pairsForUnionForEachDownPass[downPass];
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
					std::cerr << "from line: " << __LINE__<< ": downPass = " << GetPatData.toSymbol(downPass) << " INTERSECTIONS:\n";
					const VecMaskPair & forIntersections = GetPatData.pairsForIntersectionForEachDownPass[downPass];
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
			if (downPass == GetPatData.lastBitField)
				break;
			assert(downPass < GetPatData.lastBitField);
		}
		if (scObserved)
			obsmaxparscore = currscore;
	}
	if (obsmaxparscore < maxparscore) {
		std::cerr << "from line: " << __LINE__<< '\n'; std::cerr << "maxParsScore  = " << maxparscore << " obsMaxParsScore = " << obsmaxparscore << "\n";
		this->byParsScore.resize(obsmaxparscore + 1);
	}
}

void ProbInfo::calculate(const ProbInfo & leftPI, double leftEdgeLen,
					               const ProbInfo & rightPI, double rightEdgeLen,
					               TiMatFunc fn, const MTInstance & instance) {

  fn(leftEdgeLen, GetPatData.firstMatVec.GetAlias());
  const double *** leftPMatVec = const_cast<const double ***>(GetPatData.firstMatVec.GetAlias());
  fn(rightEdgeLen, GetPatData.secMatVec.GetAlias());
  const double *** rightPMatVec = const_cast<const double ***>(GetPatData.secMatVec.GetAlias());

  const unsigned int leftMaxP = leftPI.getMaxParsScore();
  const unsigned rightMaxP = rightPI.getMaxParsScore();
  const unsigned maxParsScore = 1 + leftMaxP + rightMaxP;
  unsigned int nRat = GetPatData.GetNumRates();
  unsigned int nStat = GetPatData.GetNumStates();

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
    for (std::vector<BitField>::const_iterator scIt = GetPatData.singleStateCodes.begin();
         scIt != GetPatData.singleStateCodes.end();
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
      const unsigned numStatesInMask = GetPatData.getnstates(downPass);
      if (numStatesInMask - 1 <= currScore) { // we cannot demand 3 states seen, but only 1 parsimony change... (all the probs will be zero, so we can skip them)

        MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];

        if (GetPatData.getnstates(downPass) > 1) {
          std::cerr << "from line: " << __LINE__<< ": downPass = " << (int)downPass << " " << GetPatData.toSymbol(downPass) << " UNIONS:\n";
          const VecMaskPair & forUnions = GetPatData.pairsForUnionForEachDownPass[downPass];
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
          std::cerr << "from line: " << __LINE__<< ": downPass = " << GetPatData.toSymbol(downPass) << " INTERSECTIONS:\n";
          const VecMaskPair & forIntersections = GetPatData.pairsForIntersectionForEachDownPass[downPass];
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
      if (downPass == GetPatData.lastBitField)
        break;
      assert(downPass < GetPatData.lastBitField);
    }
    if (scObserved)
      obsMaxParsScore = currScore;
  }
  if (obsMaxParsScore < maxParsScore) {
    std::cerr << "from line: " << __LINE__<< '\n'; std::cerr << "maxParsScore  = " << maxParsScore << " obsMaxParsScore = " << obsMaxParsScore << "\n";
    this->byParsScore.resize(obsMaxParsScore + 1);
  }

}

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
		std::vector<double>::const_iterator wtIt = GetPatData.categStateProb.begin();
		std::vector<double>::const_iterator pIt = pVec->begin();
		for (; wtIt != GetPatData.categStateProb.end(); ++wtIt, ++pIt) {
			assert(pIt != pVec->end());
			patClassProb += (*wtIt) * (*pIt);
		}
		for (std::vector<BitField>::const_iterator scIt = GetPatData.singleStateCodes.begin();
			scIt != GetPatData.singleStateCodes.end(); ++scIt) {
			const BitField downPass = *scIt;
			this->probsByStepsThenObsStates[0][downPass] = patClassProb;
			}

		//Other Patterns
		for (unsigned i = 1; i <= maxNumSteps; ++i) {
			const ProbForParsScore & fps = rootProbInfo.getByParsScore(i);
			for (BitField downPass = 1; ; ++downPass) {
				for (BitField obsStates = 1; ; ++obsStates) {
					double patClassProb = 0.0;
					std::vector<double>::const_iterator wtIt = GetPatData.categStateProb.begin();
					const std::vector<double> * pVec = fps.getProbsForDownPassAndObsMask(downPass, obsStates);
					if (pVec) {
						std::vector<double>::const_iterator pIt = pVec->begin();
						for (; wtIt != GetPatData.categStateProb.end(); ++wtIt, ++pIt) {
							assert(pIt != pVec->end());
							patClassProb += (*wtIt) * (*pIt);
						}
					}
					this->probsByStepsThenObsStates[i][obsStates] += patClassProb;
					if (obsStates == GetPatData.lastBitField)
						break;
				}
				if (downPass == GetPatData.lastBitField)
					break;
			}
		}
	} else {
		const unsigned maxNumSteps = rootProbInfo.getMaxParsScore();
		std::vector<double> emptyRow(GetPatData.lastBitField + 1.0, 0.0);
		this->probsByStepsThenObsStates.resize(maxNumSteps + 1, emptyRow);
		const ProbForParsScore & constFPS = rootProbInfo.getByParsScore(0);
		for (std::vector<BitField>::const_iterator scIt = GetPatData.singleStateCodes.begin();
			scIt != GetPatData.singleStateCodes.end();
			++scIt) {
			const BitField downPass = *scIt;
			const std::vector<double> * pVec = constFPS.getProbsForDownPassAndObsMask(downPass, downPass);
			assert(pVec);
			double patClassProb = 0.0;
			std::vector<double>::const_iterator wtIt = GetPatData.categStateProb.begin();
			std::vector<double>::const_iterator pIt = pVec->begin();
			for (; wtIt != GetPatData.categStateProb.end(); ++wtIt, ++pIt) {
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
						std::vector<double>::const_iterator wtIt = GetPatData.categStateProb.begin();
						const std::vector<double> * pVec = fps.getProbsForDownPassAndObsMask(downPass, obsStates);
						if (pVec) {
							std::vector<double>::const_iterator pIt = pVec->begin();
							for (; wtIt != GetPatData.categStateProb.end(); ++wtIt, ++pIt) {
								assert(pIt != pVec->end());
								patClassProb += (*wtIt) * (*pIt);
							}
						}
						this->probsByStepsThenObsStates[i][obsStates] += patClassProb;
						if (obsStates == GetPatData.lastBitField)
							break;
					}
					if (downPass == GetPatData.lastBitField)
						break;
				}
			}
	}
}


// Traverse the tree (postorder) and calculate pattern class probabilities
//
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
          if (GetPatData.isMkvSymm) {
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
            if (GetPatData.isMkvSymm) {
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

// split into funcs for different uninformative patterns
void classifyData(MTInstance & instance, const BitFieldMatrix & bmat, const int * pwPtr, PatternSummary *summ)
{
	Node * nd = instance.tree.GetRoot();
	PostorderForNodeIterator poTrav = postorder(nd);
	Arc arc = poTrav.get();
	NodeIDToParsInfo nodeIDToParsInfo;
	const ParsInfo * rootParsInfo = 0L;
	assert(GetPatData.zeroVec.size() == bmat[0].size());
	assert(arc.toNode);
	while(arc.toNode) {
		std::vector<Node *> children = arc.fromNode->GetChildren();
		const unsigned numChildren = children.size();
		NodeID currNdID(arc.fromNode, 0);
		ParsInfo & currParsInfo = nodeIDToParsInfo[currNdID];
		if (numChildren == 0) {
			const unsigned taxInd = arc.fromNode->GetNumber();
			currParsInfo.calculateForTip(bmat.at(taxInd), instance);
		} else {
			Node * leftNd = children[0];
			Node * rightNd = children[1];
			NodeIDToParsInfo::iterator leftPIIt = nodeIDToParsInfo.find(NodeID(leftNd, 0));
			NodeIDToParsInfo::iterator rightPIIt = nodeIDToParsInfo.find(NodeID(rightNd, 0));
			currParsInfo.calculateForInternal(leftPIIt->second, rightPIIt->second);

			if (numChildren > 2) {
				NodeID lastNdID(nd, 1);
				ParsInfo & lastParsInfo = nodeIDToParsInfo[currNdID];
				Node * lastNd = children.at(2);
				NodeIDToParsInfo::iterator lastPIIt= nodeIDToParsInfo.find(NodeID(lastNd, 0));
				lastParsInfo.calculateForInternal(currParsInfo, lastPIIt->second);
				rootParsInfo = &currParsInfo;
			}
			else if (arc.toNode == nd)
				rootParsInfo = &currParsInfo;
		}
	}
	assert(rootParsInfo != 0L);
	if (summ) {
		summ->clear();
		for (unsigned p = 0; p < rootParsInfo->size(); ++p) {
			unsigned toAdd = (pwPtr != 0L ? (unsigned)pwPtr[p] : 1);
			summ->incrementCount(rootParsInfo->score[p], rootParsInfo->allSeen[p], toAdd);
		}
	}
}

} // namespace
