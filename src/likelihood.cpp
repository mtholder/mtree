#include "mt_tree.h"
#include <algorithm>
using namespace std;
namespace mt {

void doAnalysis(Tree &tree, CharModel &cm)
{
    NodeIterator *pnit = postorder(tree.GetRoot());
    try {
        Node * c = pnit->get();
        assert(c);
        unsigned partIndex = 0;
        while (c) {
            std::cout << c->number<< "\n";
            const double edgeLen = c->GetEdgeLen();
            if (c->IsLeaf()) {
                const LeafCharacterVector * data = (const LeafCharacterVector *) c->GetData(partIndex);
                LeafWork * work = (LeafWork *) c->GetWork(partIndex);
                cm.fillLeafWork(data, work, edgeLen);
            } else {
                Node * left = c->leftChild;
                Node * n = left->rightSib;
                const double * c1 = left->GetCLA(partIndex);
                const double * c2 = n->GetCLA(partIndex);
                InternalNodeWork * work = (InternalNodeWork *) c->GetWork(partIndex);
                cm.fillInternalWork(c1, c2, work, edgeLen);
            }
            c = pnit->next();
        }
    } catch(...) {
        delete pnit;
        throw;
    }
    delete pnit;

}


void CharModel::fillLeafWork(const LeafCharacterVector *data, LeafWork *work, double edgeLen) {
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
void CharModel::fillInternalWork(const double * cla1, const double *cla2, InternalNodeWork *, double edgeLen) {
//    const double * this->calcTransitionProb(edgeLen);
}
}

