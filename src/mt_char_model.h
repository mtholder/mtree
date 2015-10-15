#if !defined(__CHAR_MODEL_H__)
#define __CHAR_MODEL_H__
#include <vector>
#include <cmath>
#include <string>
#include <utility>
#include "ncl/nxsallocatematrix.h"
#include "pattern_class.h"
#include "mt_data.h"

namespace mt {

class LeafCharacterVector;


// stores alphas, gamma rates, numstates, and other params specfic to each partition
// for numerical optimization
class ModelParams{
    public:
        ModelParams(unsigned numStates, unsigned modelType)
            :nStates(numStates),
            model(modelType) {
            }

        void initializeModel(unsigned nsts, unsigned mdl) {
          createGammas(0.5, nsts);
        }
        void createGammas(double a, int cats);
        double GetAlpha() const {
          return this->alpha;
        }
        unsigned getnsts() {
          return this->nStates;
        }
        std::vector<double> GetGammas(){
          return this->gammaRates;
        }
        double getGammaRate(std::size_t pos) {
          return this->gammaRates[pos];
        }
    private:
        unsigned nStates;
        unsigned model;
        double alpha;
        std::vector<double> gammaRates;
};


class CharModel {
    public:
        CharModel(unsigned maxNumStates,
                  unsigned numRateCats, 
                  //unsigned ascBiasType,
                  unsigned numparts)
            :maxNStates(maxNumStates),
            nRateCats(numRateCats),
            //ascType(ascBiasType),
            rates(1, 1.0),
            rateProb(1, 1.0) {
        }
        virtual ~CharModel() {
        }
        virtual const double * GetRateCatProb()const {
            return &(rateProb[0]);
        }
        virtual unsigned GetNumRates() const {
            return nRateCats;
        }
        virtual unsigned GetNumStates() const {
            return maxNStates;
        }

        virtual double GetRate(int pos) const {
            return rates[pos];
        }
        ModelParams & GetModelParams(int model) {
            return modelList[model];
        }
        void createGammas(double a, int cats);

        // Members for PhyPatClassProb calculations
        bool isMkvSymm;
        ScopedDblThreeDMatrix firstMatVec;
        ScopedDblThreeDMatrix secMatVec;
        unsigned pVecLen;
        std::string alphabet;
        std::vector<unsigned> zeroVec;
        BitField lastBitField;
        std::vector<double> categStateProb;
        std::vector<BitField> singleStateCodes;
        std::vector<BitField> stateIndexToStateCode;
        std::vector<unsigned> stateCodeToNumStates;
        std::vector<BitField> multiStateCodes;
        std::map<BitField, std::string> stateCodesToSymbols;
        VMaskToVecMaskPair pairsForIntersectionForEachDownPass;
        VMaskToVecMaskPair pairsForUnionForEachDownPass;
        BitFieldMatrix statesSupersets;
        unsigned getnstates(BitField mask) const {
          return stateCodeToNumStates[mask];
        }
        const std::string & toSymbol(BitField sc) const {
			    return this->stateCodesToSymbols.find(sc)->second;
		    }
        stateSetContainer::const_iterator stateSetBegin() const  {
            return possObsStateSet.begin();
        }
        stateSetContainer::const_iterator stateSetEnd() const  {
            return possObsStateSet.end();
        }

        virtual void alterRateFreq(unsigned position, double value){
          this->rates[position] = value;
        }


        virtual const double * GetRootStateFreq() const = 0;
        virtual double sumLnL(const double * cla, const double * patternWt, unsigned numChars) const;
        virtual void fillLeafWork(const LeafCharacterVector *, double * claElements, double * cla, double edgeLen, unsigned numChars);
        virtual double * calcTransitionProb(double edgeLen) = 0;
        virtual void conditionOnSingleEdge(const double *beforeEdge, double * afterEdge, double edgeLen, unsigned numChars);
        virtual void initModels(unsigned numParts, unsigned modelType);
        virtual void optimizeParams(unsigned modelType){
        }
    protected:
        unsigned maxNStates;
        unsigned nRateCats;
        std::vector<double> rates;
        std::vector<double> rateProb;
        unsigned nParts;
        //unsigned ascType;
        std::vector<double> patternWeights;
        std::vector<ModelParams> modelList;
        stateSetContainer possObsStateSet;
};

class MkCharModel: public CharModel {
    public:
        MkCharModel(unsigned numStates, unsigned numRateCats, unsigned numParts)
            :CharModel(numStates, numRateCats, numParts),
            rootStateFreq(numStates, 1.0/double(numStates)) {
        }
        virtual ~MkCharModel() {
        }
        virtual const double * GetRootStateFreq() const {
            return &(rootStateFreq[0]);
        }
        //virtual double sumLnL(const double * cla, const double * patternWt, unsigned numChars ) const;
        virtual double * calcTransitionProb(double edgeLen) {
          for(auto mi = 0U; mi < nParts; mi++) {
            unsigned nsts = modelList[mi].getnsts();
            const double fns = double(nsts);
            const double fnsmo = fns - 1.0;
            const unsigned nsSq = nsts*nsts;
            for (auto ri = 0U; ri < nRateCats; ++ri) {
                const double eb = modelList[mi].getGammaRate(ri)*edgeLen;
                const double e = std::exp(-fns*eb/fnsmo);
                const double diffP = 1.0/fns - e/fns;
                const double noDiffP = 1.0/fns +  fnsmo * e/fns;
                for (auto fi = 0U; fi < maxNStates; ++fi) {
                    for (auto ti = 0U; ti < maxNStates; ++ti) {
                        const double p = (fi == ti ? noDiffP : diffP);
                        probMat[mi*ri*nsSq + ri*nsSq + maxNStates*fi + ti] = p;
                    }
                }
            }
          }
          return &(probMat[0]);
        }
    private:
        std::vector<double> probMat;
        std::vector<double> rootStateFreq;
};


//Implementation of the general covarion model is taken from Galtier (2001)
//using the approximation for transition probabilities in equations 8,9,and 10
class MkCovarCharModel: public MkCharModel {
  public:
      MkCovarCharModel(unsigned numStates, unsigned numRateCats, unsigned numParts)
          :MkCharModel(numStates, numRateCats, numParts),
          probMat(numStates*numStates*numRateCats*numRateCats, 0.0),
          rootStateFreq(numStates, 1.0/double(numStates)) {
          }
          virtual ~MkCovarCharModel() {
          }
          virtual const double * GetRootStateFreq() const {
            return &(rootStateFreq[0]);
          }

          virtual double * calcTransitionProb(double edgeLen) {

            const double fns = double(nStates);
            const double fnsmo = fns - 1.0;
            const unsigned nsSq = nStates*nStates;

            for (auto ri = 0U; ri < nRateCats; ++ri) {
              for (auto rj = 0U; rj < nRateCats; ++rj) {
                const double eb = rates[ri]*edgeLen;
                const double eb1 = switchRate*edgeLen;
                double var1, var2, var3, probRateChange;
                if (ri == rj) {
                  var1 = 2.0*(1.0 - rates[ri]);
                  var2 = 1.0 + (double(nRateCats) - 2.0) - rates[ri];
                  var3 = double(nRateCats);
                  probRateChange = std::exp(-eb1) + (1.0 - std::exp(-eb1))/double(nRateCats);

                } else {
                    var1 = 2.0 - rates[ri] - rates[rj];
                    var2 = 1.0 - rates[ri] - rates[rj];
                    var3 = 0;
                    probRateChange = (1.0 - std::exp(-eb1))/double(nRateCats);
                }

                  const double numerator = (var1 + var2*eb)*std::exp(-eb) + eb - var1;
                  const double denominator = eb*(1.0 + (var3 - 1.0)*std::exp(-eb));
                  const double avgRate = numerator/denominator;

                  const double ebavg = avgRate*edgeLen;
                  const double e = std::exp(-fns*ebavg/fnsmo);
                  const double diffP = 1.0/fns - e/fns;
                  const double noDiffP = 1.0/fns +  fnsmo * e/fns;

                for (auto fi = 0U; fi < nStates; ++fi) {
                    for (auto ti = 0U; ti < nStates; ++ti) {
                        const double p = (fi == ti ? noDiffP : diffP);
                        probMat[ri*nRateCats*nsSq + rj*nsSq+ nStates*fi + ti] = p*probRateChange;
                      }
                    }
                }
              }
              return &(probMat[0]);
            }

  private:
      std::vector<double> probMat;
      std::vector<double> rootStateFreq;
      double switchRate;
};

class MkVarNoMissingAscCharModel: public MkCharModel {
    public:
        MkVarNoMissingAscCharModel(unsigned numStates, unsigned numRateCats)
            :MkCharModel(numStates, numRateCats) {
        }
        virtual ~MkVarNoMissingAscCharModel() {
        }
        virtual double sumLnL(const double * cla, const double * patternWt, unsigned numChars ) const;
};

class MKVarMissingAscCharModel: public MkCharModel {  public:
      MkVarMissingAscCharModel(unsigned numStates, unsigned numRateCats)
          :MkCharModel(numStates, numRateCats) {
      }
      virtual ~MkVarMissingAscCharModel() {
      }
      virtual double sumLnL(const double * cla, const double * patternWt, unsigned numChars ) const;
};

class MkVarParsInfNoMissingModel: public MkCharModel {
    public:
        MkVarParsInfNoMissingModel(unsigned numStates, unsigned numRateCats)
            :MkCharModel(numStates, numRateCats) {
        }
        virtual ~MkVarParsInfNoMissingModel() {
        }
        virtual double sumLnL(const double * cla, const double * patternWt, unsigned numChars) const;
};

// Character model for case with no missing data, and only parsimony-informative patterns
class MkVarParsInfMissingModel: public MkCharModel {
    public:
        MkVarParsInfMissingModel(unsigned numStates, unsigned numRateCats)
            :MkCharModel(numStates, numRateCats) {
        }
        virtual ~MkVarParsInfMissingModel() {
        }
        virtual double sumLnL(const double * cla, const double * patternWt, unsigned numChars) const;
};


} //namespace
#endif
