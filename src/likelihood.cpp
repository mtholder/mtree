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
void doAnalysis(Tree &tree, CharModel &cm)
{
    PostorderForNodeIterator pnit = postorder(tree.GetRoot());
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

}

inline void assign(double value, double * dest, unsigned n) {
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
            for (auto ri = 0; ri < nRateCats; ++ri) {
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

void CharModel::conditionOnSingleEdge(const double * cla1, double * afterEdge, double edgeLen, unsigned numChars) {
    //const double * tiprob = this->calcTransitionProb(edgeLen);
    
}
/*
void CharModel::fillInternalWork(const double * cla1, const double *cla2, InternalNodeWork *, double edgeLen) {
//    const double * this->calcTransitionProb(edgeLen);
}
*/
}

