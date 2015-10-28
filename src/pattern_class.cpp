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

// TEMP - throws exception for missing data
// return a string of symbols for each state with length = num states
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
    //NxsCharactersBlock::DataTypesEnum inDatatype = mapper.GetDatatype();
    const unsigned nStates =  mapper.GetNumStates();
    //assert(nStates <= MAX_NUM_STATES);
    for (NxsDiscreteStateCell i = 0; i < (NxsDiscreteStateCell)nStates; ++i)
        assert(mapper.GetStateSetForCode(i).size() == 1);
    const std::string fundamentalSymbols = mapper.GetSymbols();
    assert(fundamentalSymbols.length() == nStates);
    const unsigned nTaxa = charsBlock.GetNTax();
    const auto includedNChar = toInclude->size();
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



// Class to store node data for uninformative case
// NodeDataStructure in Uninformative_case.cpp
class NodeInfo {
    public:
        NodeInfo(unsigned numStates)
          :numLeaves(0),
          isMissing(false) {
          int len = 1 << numStates; // bitwise shift because character states as bitfields
          ProbForObsStateSet dummy (numStates);
          probVec.assign(len, dummy);
        }

        ProbForObsStateSet & getForObsStateSet(int obs) {
            return probVec.at(obs);
        }

        void copyProbVec(NodeInfo *info) {
          for(auto i = 0U; i < probVec.size(); i++) {
            probVec[i] = info->getForObsStateSet(i);
          }
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

int countBits(int x);
int convertIndexToBit(int ind);
int convertBitToIndex(int i);
std::vector<int> subsetsContainingGivenState(int fullSet, int givenState);
std::vector<int> subsetsOfGivenSize(int obsStSet, int numBits);
int getNextCommStSet(const int obsStSet, int i);
double pclassCalcTransitionProb(int ancIndex, int i, double edgeLen, MTInstance & instance);
double calcProbOfSubtreeForObsStSetAndComm(NodeInfo * subtreeInfo, int ancIndex, int obsBits, int commonStates, double edgeLen, MTInstance &instance);
double calcProbOfSubtreeForObsStSetNoRepeated(NodeInfo * subtreeInfo, int ancIndex, int obsBits, double edgeLen, MTInstance &instance);
void calcUninformativePatterns(MTInstance & instance);

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

double pclassCalcTransitionProb(int ancIndex, int i, double edgeLen, MTInstance & instance, unsigned model){
  double * tiVec = GetPatData(model).calcTransitionProb(edgeLen);
  int nStates = GetPatData(model).GetNumStates();
  return tiVec[ancIndex*nStates + i];
}

double calcProbOfSubtreeForObsStSetAndComm(NodeInfo * subtreeInfo, int ancIndex, int obsBits, int commonStates, double edgeLen,
                                           MTInstance &instance, unsigned model) {
  double p = 0.0;
  ProbForObsStateSet & childProbSet = subtreeInfo->getForObsStateSet(obsBits);
  std::vector<double> & childProb = childProbSet.getProbForCommState(commonStates);
  for(auto i = 0U; i < GetPatData(model).GetNumStates(); i++)  {
    double transProb = pclassCalcTransitionProb(ancIndex, i, edgeLen, instance, model);
    double partialLike = childProb[i];
    double x = transProb * partialLike;
    p += x;
  }
  return p;
}

double calcProbOfSubtreeForObsStSetNoRepeated(NodeInfo * subtreeInfo, int ancIndex, int obsBits, double edgeLen, MTInstance &instance, unsigned model){
  return calcProbOfSubtreeForObsStSetAndComm(subtreeInfo, ancIndex, obsBits, -1, edgeLen, instance, model);
}

#if 0
// Traverse the tree (postorder) and calculate pattern class probabilities
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

#endif

void cleanVirtualEdgeLens(Node * root) {
  PostorderForNodeIterator poTrav = postorder(root);
  Arc arc = poTrav.get();
  do {
    arc.fromNode->SetVEdgeLen(0.0);
    arc = poTrav.next();
  } while(arc.toNode);
}

// calculate probabilities of uninformative patterns
// for one partition, for subtree rooted at nd
NodeInfo * calcUninformativePatterns(MTInstance & instance, Node * nd, unsigned charIndex, unsigned model)
{
  //Node * nd = instance.tree.GetRoot();
  PostorderForNodeIterator poTrav = postorder(nd);
  Arc arc = poTrav.get();
  std::map<Node *, NodeInfo *> nodeToInfoMap;
  unsigned numStates = GetPatData(model).GetNumStates();
  NodeInfo * currNdInfo = 0L;
  assert(arc.toNode);
  while(arc.toNode) {
    Node * currNd = arc.fromNode;
    std::vector<Node *> children = currNd->GetChildren();
    const auto numChildren = children.size();
    currNdInfo = new NodeInfo(numStates); // should this be on heap or stack?
    NodeID currNdID(currNd, 0);
    nodeToInfoMap[currNd] = currNdInfo;
    if (numChildren == 0) {
      unsigned li = currNd->GetNumber();
      if (instance.partMat.GetLeafCharacters(model,li)->GetCharVec()[charIndex] == numStates) {
        currNdInfo->setIsMissing(true);
      }
      for(auto i = 0U; i < numStates; i++) {
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
      // if (numChildren == 2)
      Node * leftChild = children[0];
      NodeInfo * leftNdInfo = nodeToInfoMap[leftChild];
      Node * rightChild = children[1];
      NodeInfo * rightNdInfo = nodeToInfoMap[rightChild];
      if(leftNdInfo->getIsMissing() || rightNdInfo->getIsMissing()) {
        if(leftNdInfo->getIsMissing() && rightNdInfo->getIsMissing()) {
          currNdInfo->setIsMissing(true);
          arc = poTrav.next();
          continue;
        } else {
          // add edge lengths and progress to next node
          if(leftNdInfo->getIsMissing()) {
            currNd->SetVEdgeLen(rightChild->GetEdgeLen());
            currNdInfo->setNumLeaves(rightNdInfo->getNumLeaves());
            currNdInfo->copyProbVec(rightNdInfo);
            arc = poTrav.next();
            continue;
          } else { // right is missing
            currNd->SetVEdgeLen(leftChild->GetEdgeLen());
            currNdInfo->setNumLeaves(leftNdInfo->getNumLeaves());
            currNdInfo->copyProbVec(leftNdInfo);
            arc = poTrav.next();
            continue;
          }
        }
      } else {
      currNdInfo->setNumLeaves(leftNdInfo->getNumLeaves() + rightNdInfo->getNumLeaves());

      stateSetContainer::const_iterator ssCit = GetPatData(model).stateSetBegin();
      for (; ssCit != GetPatData(model).stateSetEnd(); ssCit++) {
        const int & obsStSet = *ssCit;
        int common = -1;   // this is 111... in bits
        int numObsSt = countBits(obsStSet);

        while(common>-2) {
          ProbForObsStateSet & currNdProbSet = currNdInfo->getForObsStateSet(obsStSet);
          std::vector<double> & currNdProbVec = currNdProbSet.getProbForCommState(common);

          if(common == -1) {
            if (currNdInfo->getNumLeaves() == numObsSt) {
              for (auto anc = 0U; anc < numStates; anc++) {
                std::cerr << "ObsStSet " << obsStSet << '\n';
                currNdProbVec[anc] = 0.0;
                std::vector<int> leftObsStSets = subsetsOfGivenSize(obsStSet, leftNdInfo->getNumLeaves());
                for (auto j = 0U; j < leftObsStSets.size(); j++) {
                  int leftObsStSet = leftObsStSets[j];
                  int rightObsStSet = obsStSet - leftObsStSet;

                  double leftProb, rightProb;
                  double leftedgeLen = leftChild->GetEdgeLen() + leftChild->GetVEdgeLen();
                  if(leftNdInfo->getNumLeaves() == 1) {
                    leftProb = pclassCalcTransitionProb(anc, convertBitToIndex(leftObsStSet), leftedgeLen, instance, model);
                  } else {
                    leftProb = calcProbOfSubtreeForObsStSetNoRepeated(leftNdInfo, anc, leftObsStSet, leftedgeLen, instance, model);
                  }

                  double rightEdgeLen = rightChild->GetEdgeLen() + rightChild->GetVEdgeLen();
                  if(rightNdInfo->getNumLeaves() == 1) {
                    rightProb = pclassCalcTransitionProb(anc, convertBitToIndex(rightObsStSet), rightEdgeLen, instance, model);
                  } else {
                    rightProb = calcProbOfSubtreeForObsStSetNoRepeated(rightNdInfo, anc, rightObsStSet, rightEdgeLen, instance, model);
                  }
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;
                  }
                }
              }
            } else {

              int commonBits = convertIndexToBit(common);
              for (auto anc = 0U; anc < numStates; anc++) {

                currNdProbVec[anc] = 0.0;
                int leftCommSt, rightCommSt;
                leftCommSt = common;
                rightCommSt = common;
                std::vector<int> obsStSetsWithComm = subsetsContainingGivenState(obsStSet, commonBits);

                for(auto j=0U; j < obsStSetsWithComm.size(); j++) {
                  int leftObsStSet = obsStSetsWithComm[j];
                  int rightObsStSet = obsStSet - leftObsStSet + commonBits;

                  double leftProb, rightProb;
                  double leftEdgeLen = leftChild->GetEdgeLen() + leftChild->GetVEdgeLen();
                  leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, common, leftEdgeLen, instance, model);
                  double rightEdgeLen = rightChild->GetEdgeLen() + rightChild->GetVEdgeLen();
                  rightEdgeLen = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, common, rightEdgeLen, instance, model);
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;
                }

                leftCommSt = -1;
                rightCommSt = common;
                //add probability when only right common, left not repeated
                for(auto j = 0U; j < obsStSetsWithComm.size(); j++) {
                  int rightObsStSet = obsStSetsWithComm[j];
                  int leftObsStSet = obsStSet - rightObsStSet + commonBits;
                  double leftProb, rightProb;
                  double leftEdgeLen = leftChild->GetEdgeLen() + leftChild->GetVEdgeLen();
                  leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, -1, leftEdgeLen, instance, model);
                  double rightEdgeLen = rightChild->GetEdgeLen() + rightChild->GetVEdgeLen();
                  rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, common, rightEdgeLen, instance, model);
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;

                  //Now consider when the left is not displayed by commonBits as observed States
                  leftObsStSet = obsStSet - rightObsStSet;
                  if(leftObsStSet != 0) {
                    leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, -1, leftEdgeLen, instance, model);
                    jointNdProb = leftProb * rightProb;
                    currNdProbVec[anc] += jointNdProb;
                    }
                }

                leftCommSt = common;
                rightCommSt = -1;
                //add probability when only left common
                for(auto j = 0U; j < obsStSetsWithComm.size(); j++) {
                  int rightObsStSet = obsStSetsWithComm[j];
                  int leftObsStSet = obsStSet - rightObsStSet + commonBits;
                  double leftProb, rightProb;
                  double leftEdgeLen = leftChild->GetEdgeLen() + leftChild->GetVEdgeLen();
                  leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, common, leftEdgeLen, instance, model);
                  double rightEdgeLen = rightChild->GetEdgeLen() + rightChild->GetVEdgeLen();
                  rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, -1, rightEdgeLen, instance, model);
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;

                  //Now consider when the right is not displayed by commonBits as observed States
                  rightObsStSet = obsStSet - leftObsStSet;
                  if(rightObsStSet != 0) {
                    rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, -1, rightEdgeLen, instance, model);
                    jointNdProb = leftProb * rightProb;
                    currNdProbVec[anc] += jointNdProb;
                  }
                }

                leftCommSt = -1;
                rightCommSt = -1;
                //add probability when neither common
                for(auto j = 0U; j < obsStSetsWithComm.size(); j++) {
                  int rightObsStSet = obsStSetsWithComm[j];
                  int leftObsStSet = obsStSet - rightObsStSet + commonBits;
                  double leftProb, rightProb;
                  double leftEdgeLen = leftChild->GetEdgeLen() + leftChild->GetVEdgeLen();
                  leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNdInfo, anc, leftObsStSet, -1, leftEdgeLen, instance, model);
                  double rightEdgeLen = rightChild->GetEdgeLen() + rightChild->GetVEdgeLen();
                  rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNdInfo, anc, rightObsStSet, -1, rightEdgeLen, instance, model);
                  double jointNdProb = leftProb * rightProb;
                  currNdProbVec[anc] += jointNdProb;
                }
              }
            }
            common = getNextCommStSet(obsStSet, common);
          }
        }
      arc = poTrav.next();
    }
    }
  }
    return currNdInfo;
    //do something else
} // calcUninformativePatterns function end


} // namespace
