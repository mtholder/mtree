#include "mt_instance.h"
#include "mt_data.h"
#include "mt_tree.h"
#include "mt_tree_traversal.h"
#include "parsimony.h"
#include "pattern_class.h"

#include <cassert>
#include <vector>

namespace mt {
// Functions for parsimony calcs used in likelihood search and pattern class calcs
// Using Fitch's algorithm
unsigned int calcParsScore(MTInstance & instance);
#if 0
// Set downpass state to taxon data,
void ParsInfo::calculateForTip(const BitFieldRow & data, MTInstance & instance) {
    assert(instance.GetCharModel(0).zeroVec.size() >= data.size());
    this->numPatterns = data.size();
    this->downPass = &data[0];
    this->allSeen = &data[0];
    this->score = &instance.GetCharModel(0).zeroVec[0];
}
void ParsInfo::calculateForInternal(ParsInfo & leftData, ParsInfo & rightData) {
    this->numPatterns = leftData.size();
    assert(numPatterns == rightData.size());
    this->downPassOwned.resize(numPatterns);
    this->downPass = &this->downPassOwned[0];
    this->allSeenOwned.resize(numPatterns);
    this->allSeen = &this->allSeenOwned[0];
    this->scoreOwned.resize(numPatterns);
    this->score = &this->scoreOwned[0];
    const BitField * leftDown = leftData.downPass;
    const BitField * leftAll = leftData.allSeen;
    const std::size_t * leftScore = leftData.score;
    const BitField * rightDown = rightData.downPass;
    const BitField * rightAll = rightData.allSeen;
    const std::size_t * rightScore = rightData.score;
    for (unsigned i = 0; i < numPatterns; ++i, ++leftDown, ++leftAll, ++leftScore, ++rightDown, ++rightAll, ++rightScore) {
        const BitField lrIntersection = (*leftDown) & (*rightDown);
        const std::size_t accumulScore = *leftScore + *rightScore;
        this->allSeenOwned[i] = (*leftAll) | (*rightAll);
        if (lrIntersection == BitField(0)) {
            const BitField lrUnion = (*leftDown) | (*rightDown);
            this->downPassOwned[i] = lrUnion;
            this->scoreOwned[i] = 1 + accumulScore;
        }
        else {
            this->downPassOwned[i] = lrIntersection;
            this->scoreOwned[i] = accumulScore;
        }
    }
}

void ParsInfo::write(std::ostream & o) const {
    unsigned totalScore = 0;
    for (unsigned i = 0; i < numPatterns; ++i) {
        o << "pattern = " << i << " score = " << (int)score[i] << " downpass = " << (int)downPass[i] << " states = " << (int)allSeen[i] << '\n';
        totalScore += score[i];
    }
    o << "totalScore = " << totalScore << '\n';
}

// Return parsimony score of tree
unsigned int calcParsScore(MTInstance & instance) {
    unsigned int score = 0;
    Node * nd = instance.tree.GetRoot();
    NodeID currNdId(nd, 0);
    Node * vRoot = nd->leftChild->rightSib;
    PostorderForNodeIterator postorderTrav = postorder(vRoot);
    Arc ptArc = postorderTrav.get();
    while (ptArc.toNode) {
        //score += ...
        ptArc = postorderTrav.next();
    }
    return score;
}
#endif
} //namespace
