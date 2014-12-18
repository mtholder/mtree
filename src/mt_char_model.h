#if !defined(__CHAR_MODEL_H__)
#define __CHAR_MODEL_H__
#include <vector>
#include <cmath>
namespace mt {

class LeafCharacterVector;

class CharModel {
    public:
        CharModel(unsigned numStates, unsigned numRateCats)
            :nStates(numStates),
            nRateCats(numRateCats),
            rates(1, 1.0), 
            rateProb(1, 1.0) {
        }
        virtual ~CharModel() {
        }
        virtual const double * GetRateCatProb() const {
            return &(rateProb[0]);
        }
        virtual unsigned GetNumRates() const {
            return nRateCats;
        }
        virtual unsigned GetNumStates() const {
            return nStates;
        }
        
        virtual const double * GetRootStateFreq() const = 0;
        virtual double sumLnL(const double * cla, const double * patternWt, unsigned numChars) const;
        virtual void fillLeafWork(const LeafCharacterVector *, double * claElements, double * cla, double edgeLen, unsigned numChars);
        virtual double * calcTransitionProb(double edgeLen) = 0;
        virtual void conditionOnSingleEdge(const double *beforeEdge, double * afterEdge, double edgeLen, unsigned numChars);
    protected:
        unsigned nStates;
        unsigned nRateCats;
        std::vector<double> rates;
        std::vector<double> rateProb;
};

class MkCharModel: public CharModel {
    public:
        MkCharModel(unsigned numStates, unsigned numRateCats)
            :CharModel(numStates, numRateCats),
            probMat(numStates*numStates*numRateCats, 0.0),
            rootStateFreq(numStates, 1.0/float(numStates)) {
        }
        virtual ~MkCharModel() {
        }
        virtual const double * GetRootStateFreq() const {
            return &(rootStateFreq[0]);
        }
        //virtual double sumLnL(const double * cla, const double * patternWt, unsigned numChars ) const;
        virtual double * calcTransitionProb(double edgeLen) {
            const double fns = float(nStates);
            const double fnsmo = fns - 1.0;
            const unsigned nsSq = nStates*nStates;
            for (auto ri = 0U; ri < nRateCats; ++ri) {
                const double eb = rates[ri]*edgeLen;
                const double e = std::exp(-fns*eb/fnsmo);
                const double diffP = 1.0/fns - e/fns;
                const double noDiffP = 1.0/fns +  fnsmo * e/fns;
                for (auto fi = 0U; fi < nStates; ++fi) {
                    for (auto ti = 0U; ti < nStates; ++ti) {
                        const double p = (fi == ti ? noDiffP : diffP);
                        probMat[ri*nsSq + nStates*fi + ti] = p;
                    }
                }
            }
            return &(probMat[0]);
        }
    private:
        std::vector<double> probMat;
        std::vector<double> rootStateFreq;
};

class MkVarNoMissingAscCharModel: public MkCharModel {
    public:
        MkVarNoMissingAscCharModel(unsigned numStates, unsigned numRateCats)
            :MkCharModel(numStates, numRateCats) {
        }
        virtual ~MkVarNoMissingAscCharModel() {
        }
        //virtual double sumLnL(const double * cla, const double * patternWt, unsigned numChars ) const;
};


} //namespace
#endif
