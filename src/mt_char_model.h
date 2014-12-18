#if !defined(__CHAR_MODEL_H__)
#define __CHAR_MODEL_H__
namespace mt {

class LeafCharacterVector;

class CharModel {
    public:
        CharModel(unsigned numStates, unsigned numRateCats)
            :nStates(numStates),
            nRateCats(numRateCats) {
        }
        virtual ~CharModel() {
        }
        virtual double sumLnL(const double * cla, unsigned numChars) const = 0;
        virtual void fillLeafWork(const LeafCharacterVector *, double * claElements, double * cla, double edgeLen, unsigned numChars);
        virtual double * calcTransitionProb(double edgeLen) = 0;
        virtual void conditionOnSingleEdge(const double *beforeEdge, double * afterEdge, double edgeLen, unsigned numChars);
    protected:
        unsigned nStates;
        unsigned nRateCats;
};
class MkVarNoMissingAscCharModel: public CharModel {
    public:
        MkVarNoMissingAscCharModel(unsigned numStates, unsigned numRateCats)
            :CharModel(numStates, numRateCats) {
        }
        virtual ~MkVarNoMissingAscCharModel() {
        }
        virtual double sumLnL(const double * cla, unsigned numChars ) const {
            return 3.2;
        }
        virtual double * calcTransitionProb(double edgeLen) {
            return nullptr;
        }
};
class MkCharModel: public CharModel {
    public:
        MkCharModel(unsigned numStates, unsigned numRateCats)
            :CharModel(numStates, numRateCats) {
        }
        virtual ~MkCharModel() {
        }
        virtual double sumLnL(const double * cla, unsigned numChars ) const {
            return 3.4;
        }
        virtual double * calcTransitionProb(double edgeLen) {
            return nullptr;
        }
};


} //namespace
#endif
