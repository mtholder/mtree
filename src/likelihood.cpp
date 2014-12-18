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
    while (c.toNode) {
        std::cout << c.fromNode->number<< "\n";
        const double edgeLen = c.GetEdgeLen();
        if (c.IsFromLeaf()) {
            const LeafCharacterVector * data = c.GetFromNdData(partIndex);
            LeafWork * work = c.GetFromNdLeafWork(partIndex);
            cm.fillLeafWork(data, work, edgeLen);
        } else {
            vector<const double *> p = c.GetPrevCLAs(partIndex);
            double * beforeArc = c.GetFromNdCLA(partIndex, false);
            pruneProductStep(p, beforeArc, c.GetLenCLA(partIndex));
            double * afterArc = c.GetFromNdCLA(partIndex, true);
            cm.conditionOnSingleEdge(beforeArc, afterArc, edgeLen, c.GetNumChars(partIndex));
        }
        c = pnit.next();
    }

}


void CharModel::fillLeafWork(const LeafCharacterVector *data, LeafWork *LeafWork, double edgeLen) {
    /* fill the summed probabilities for each state code */
#if 0
    const double * tiprob = this->calcTransitionProb(edgeLen);
    const CharStateToPrimitiveInd * s2pi = data->cs2pi;
    const unsigned numStateCodes = s2pi->GetNumStateCodes();
    const unsigned lenCLAWord = nStates*nRateCats;
    double * summedLoc = &(work->summed[0]);
    const unsigned nssq = nStates * nStates;
    for (auto sci = 0U; sci < numStateCodes; ++sci) {
        assign(0.0, lenCLAWord, summedLoc);
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
    const numChars = data->charVec.size();
    double * cla = work->cla;
    for (auto ci = 0U; ci < numChars; ++ci) {
        const char_state_t sc = data->charVec[ci];
        const double * s = &(work->summed[sc]);
        copy(s, s + lenCLAWord, cla + ci*lenCLAWord);
    }
#endif

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

