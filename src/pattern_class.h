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
typedef std::pair<const Node *, unsigned int> NodeID;
typedef std::map<BitField, std::vector<double> > MaskToProbsByState;
typedef std::map<BitField, MaskToProbsByState > MaskToMaskToProbsByState;

class ProbInfo;
typedef std::map<NodeID, ProbInfo *> NodeIDToProbInfo;

typedef std::pair<BitField, BitField> MaskPair;
typedef std::vector<MaskPair> VecMaskPair;
typedef std::map<BitField, VecMaskPair> MaskToVecMaskPair;
typedef std::vector<VecMaskPair> VMaskToVecMaskPair;

const std::vector<double> * getProbsForStatesMask(const MaskToProbsByState *, const BitField sc);
std::string convertToBitFieldMatrix(const NxsCharactersBlock & charsBlock, BitFieldMatrix & bfMat);

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

class ProbForParsScore {
    public:
        const MaskToProbsByState * getMapPtrForDownPass(const BitField sc) const {
          MaskToMaskToProbsByState::const_iterator scIt = this->byDownPass.find(sc);
          if (scIt == this->byDownPass.end())
              return 0L;
          return &(scIt->second);
        }
        const std::vector<double> * getProbsForDownPassAndObsMask(const BitField downPass, const BitField mask) const {
          return getProbsForStatesMask(this->getMapPtrForDownPass(downPass), mask);
        }

    private:
    // map from Downpass BitField => map of observed state set BitField => prob vec
      MaskToMaskToProbsByState byDownPass;
      friend class ProbInfo;
};

class MTInstance;

class ExpectedPatternSummary {
  public:
      ExpectedPatternSummary(const ProbInfo &, const MTInstance &);
      void write(std::ostream, const MTInstance &) const;
  private:
      std::vector<std::vector<double> > probsByStepsThenObsStates;
};

class ProbInfo {
  public:
      void createForTip(const MTInstance &);
      void calculateSymmetric(const ProbInfo & leftPI, double leftEdgeLen,
			    const ProbInfo & rightPI, double rightEdgeLen,
			    TiMatFunc fn, const MTInstance &);
      void calculate(const ProbInfo & leftPI, double leftEdgeLen,
					const ProbInfo & rightPI, double rightEdgeLen,
					TiMatFunc fn, const MTInstance &);
      unsigned getMaxParsScore() const {
		    assert(!this->byParsScore.empty());
		    return this->byParsScore.size() - 1;
      }
      const ProbForParsScore & getByParsScore(unsigned score) const {
		    return this->byParsScore.at(score);
      }
      unsigned getNLeavesBelow() const {
		    return nLeavesBelow;
      }
  protected:
      void addToAncProbVec(
                std::vector<double> & pVec,
                const double *** leftPMatVec, const std::vector<double> * leftProbs,
                const double *** rightPMatVec, const std::vector<double> * rightProbs,
                const MTInstance & instance);
        // declaration
      void addToAncProbVecSymmetric(std::vector<double> & pVec,
                const double *** leftPMatVec, const std::vector<double> * leftProbs,
                const double *** rightPMatVec, const std::vector<double> * rightProbs,
                const std::vector<unsigned int> & rightChildStateCodeTranslation,
                const MTInstance & instance);

      bool allCalcsForAllPairs(
                MaskToProbsByState & forCurrScoreDownPass,
                const VecMaskPair & pairVec,
                const ProbInfo & leftPI,
                const double *** leftPMatVec,
                const ProbInfo & rightPI,
                const double *** rightPMatVec,
                const unsigned accumScore,
                const bool doingIntersection,
                const MTInstance & blob);
      // data
      unsigned nLeavesBelow;
      std::vector<ProbForParsScore> byParsScore;
};

} //namespace
#endif
