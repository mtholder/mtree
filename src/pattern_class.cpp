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
// https://github.com/mtholder/PhyPatClassProb/blob/master/src/Uninformative_case.cpp

namespace mt {

//debug function
void Tree::initBrLens() {
    PostorderForNodeIterator pTrav = postorder(root);
    Arc arc = pTrav.get();
    while(arc.toNode) {
      arc.fromNode->SetEdgeLen(0.1);
      arc = pTrav.next();
    }
}

void patClassInitialize(MTInstance &instance) {
  for(int m = 0; m < instance.numPartitions; m++) {
    int numstates = GetPatData(m).GetNumStates();
    int endIndex = 1<<(numstates);
    for (int i = 1; i < endIndex; i++){
      GetPatData(m).possObsStateSet.push_back(i);
    }
  }
}

// check that categ probs sum to 1
void addCategProbs(const std::vector<double> & v, int n) {
  double sum = 0.0;
  for(int i = 0; i < n; i++){
    sum += v[i];
  }
  sum -= -1.0;
  assert(fabs(sum < 1e-6));
}

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
        void write_probVec() const {
          for (std::vector<ProbForObsStateSet>::const_iterator pvIt = probVec.begin(); pvIt != probVec.end(); pvIt++) {
            const ProbForObsStateSet &pfoss = *pvIt;
            pfoss.write_pv();
          }

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
  //std::cerr << "ObsStSet = " << obsStSet << ", i = " << i << "\n";
  int ind, binRep;
  if(i == -1) {
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

// free all NodeInfo structures above root NodeInfo
void freeNodeInfo(MTInstance &instance, Node * root, std::map<Node *, NodeInfo *> nodeToInfoMap) {
  PostorderForNodeIterator poTrav = postorder(root);
  Arc arc = poTrav.get();
  while (arc.toNode) {
    Node * currNd = arc.fromNode;
    //nodeToInfoMap[currNd]->write_probVec();
    //std::cerr << "Deleting nodeinfo\n";
    delete nodeToInfoMap[currNd];
    //std::cerr << "Deleted nodeinfo\n";
    //std::cerr << "Erasing key in map\n";
    nodeToInfoMap.erase(currNd);
    //std::cerr << "Erased key in map\n";
    arc = poTrav.next();
  }
}

double pclassCalcTransitionProb(int ancIndex, int i, double edgeLen, MTInstance & instance, unsigned model){
  double * tiVec = GetPatData(model).calcTransitionProb(edgeLen);
  int nStates = GetPatData(model).GetNumStates();
  int nRateCats = GetPatData(model).GetNumRates();
  int nssq = nStates*nStates;
  const double * rateProbs = GetPatData(model).GetRateCatProb();
  // marginalize over different rates
  double sumProb = 0.0;
  for (auto ri = 0U; ri < nRateCats; ri++) {
    sumProb += tiVec[ri*nssq + ancIndex*nStates + i]*rateProbs[ri];
  }

  return sumProb;
}

double calcProbOfSubtreeForObsStSetAndComm(NodeInfo * subtreeInfo, int ancIndex, int obsBits, int commonStates, double edgeLen,
                                           MTInstance &instance, unsigned model) {
  //std::cerr << "Entering calcProbOfSubtreeForObsStSetAndComm\n";
  double p = 0.0;
  ProbForObsStateSet & childProbSet = subtreeInfo->getForObsStateSet(obsBits);
  std::vector<double> & childProb = childProbSet.getProbForCommState(commonStates);
  for(auto i = 0U; i < GetPatData(model).GetNumStates(); i++)  {
    double transProb = pclassCalcTransitionProb(ancIndex, i, edgeLen, instance, model);
    //std::cerr << "transProb = " << transProb << "\n";
    double partialLike = childProb[i];
    if(partialLike > 1.0){
      partialLike = 0.0;
    }
    //std::cerr << "partialLike = " << partialLike << "\n";
    double x = transProb * partialLike;
    p += x;
  }
  //std::cerr << "p = " << p << "\n";
  return p;
}

double calcProbOfSubtreeForObsStSetNoRepeated(NodeInfo * subtreeInfo, int ancIndex, int obsBits, double edgeLen, MTInstance &instance, unsigned model){
  //std::cerr << "Entering calcProbOfSubtreeForObsStSetNoRepeated\n";
  return calcProbOfSubtreeForObsStSetAndComm(subtreeInfo, ancIndex, obsBits, -1, edgeLen, instance, model);
}


void cleanVirtualEdgeLens(Node * root) {
  PostorderForNodeIterator poTrav = postorder(root);
  Arc arc = poTrav.get();
  do {
    arc.fromNode->SetVEdgeLen(0.0);
    arc = poTrav.next();
  } while(arc.toNode);
}

// calculate ln probabilities of uninformative patterns
// for one partition, for subtree rooted at nd
double calcUninformativePatterns(MTInstance & instance, Node * nd, unsigned charIndex, unsigned model)
{
  //Node * nd = instance.tree.GetRoot();
  //std::cout << "Entering calcUninformativePatterns\n";
  double total = 0.0;
  PostorderForNodeIterator poTrav = postorder(nd);
  Arc arc = poTrav.get();
  std::map<Node *, NodeInfo *> nodeToInfoMap;
  unsigned numStates = GetPatData(model).GetNumStates();
  //std::cerr << "Model = " << model <<", numStates = " << numStates << "\n";
  NodeInfo * currNdInfo = 0L;
  assert(arc.toNode);
  while(arc.toNode) {
    Node * currNd = arc.fromNode;
    //std::cerr << "Branch length: " << currNd->GetEdgeLen() << "\n";
    std::vector<Node *> children = currNd->GetChildren();
    const auto numChildren = children.size();
    //std::cout << "Creating new NodeInfo, Node = " << currNd->GetNumber() << "\n";
    currNdInfo = new NodeInfo(numStates);
    NodeID currNdID(currNd, 0);
    nodeToInfoMap[currNd] = currNdInfo;
    if (numChildren == 0) {
      // Node is a leaf
      unsigned li = currNd->GetNumber();
      if (instance.partMat.GetLeafCharacters(model,li)->GetCharVec()[charIndex] == numStates) {
        //std::cout << "Data missing at character " << charIndex << "at leaf "<< li <<"\n";
        currNdInfo->setIsMissing(true);
      }
      for(auto i = 0U; i < numStates; i++) {
        int ss=1 << i;
        ProbForObsStateSet & p = currNdInfo->getForObsStateSet(ss);
        std::vector<double> & v = p.getProbForCommState(-1);
        v[i] = 1.0;
        currNdInfo->setNumLeaves(1);
      }
      arc = poTrav.next();
      continue;
    } else {

    if (numChildren != 2) {
      std::cerr << "Trees must be binary \n";
      exit(1);
    }

      // Node is internal node
      Node * leftChild = children[0];
      NodeInfo * leftNdInfo = nodeToInfoMap[leftChild];
      Node * rightChild = children[1];
      NodeInfo * rightNdInfo = nodeToInfoMap[rightChild];
      if(leftNdInfo->getIsMissing() || rightNdInfo->getIsMissing()) {
        if(leftNdInfo->getIsMissing() && rightNdInfo->getIsMissing()) {
          currNdInfo->setIsMissing(true);
          arc = poTrav.next();
          continue;
          //std::cerr << "Continuing after missing data @ internal node\n";
        } else {
          // add edge lengths and progress to next node
          if(leftNdInfo->getIsMissing()) {
            currNd->SetVEdgeLen(rightChild->GetEdgeLen());
            currNdInfo->setNumLeaves(rightNdInfo->getNumLeaves());
            currNdInfo->copyProbVec(rightNdInfo);
            arc = poTrav.next();
            continue;
            //std::cerr << "Continuing after missing data @ left child\n";
          } else { // right is missing
            currNd->SetVEdgeLen(leftChild->GetEdgeLen());
            currNdInfo->setNumLeaves(leftNdInfo->getNumLeaves());
            currNdInfo->copyProbVec(leftNdInfo);
            arc = poTrav.next();
            continue;
            //std::cerr << "Continuing after missing data @ right child\n";
          }
        }
      } else {
      currNdInfo->setNumLeaves(leftNdInfo->getNumLeaves() + rightNdInfo->getNumLeaves());

      //std::cerr << "Reading stateset iterator\n";
      stateSetContainer::const_iterator ssCit = GetPatData(model).stateSetBegin();
      //std::cerr << "Read stateset iterator\n";
      for (; ssCit != GetPatData(model).stateSetEnd(); ssCit++) {
        //std::cerr << "got this far\n";
        const int & obsStSet = *ssCit;
        //std::cerr << "got this far\n";
        int common = -1;   // this is 111... in bits
        int numObsSt = countBits(obsStSet);

        while(common>-2) {
          //std::cerr << "Common " << common << "\n";
          ProbForObsStateSet & currNdProbSet = currNdInfo->getForObsStateSet(obsStSet);
          std::vector<double> & currNdProbVec = currNdProbSet.getProbForCommState(common);

          if(common == -1) {
            if (currNdInfo->getNumLeaves() == numObsSt) {
              for (auto anc = 0U; anc < numStates; anc++) {
                //std::cerr << "ObsStSet " << obsStSet << '\n';
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

                //std::cerr << "State set size = " << obsStSetsWithComm.size() << "\n";
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
                  //std::cerr << "Right common\n";
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
                  //std::cerr << "Left common\n";
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
                  //std::cerr << "Neither common\n";
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
            //std::cerr << "got here\n";
            common = getNextCommStSet(obsStSet, common);
          }
        }
      arc = poTrav.next();
    }
    }
  }
  // iterate over statesets to sum probs for this partition + character
  stateSetContainer::const_iterator ssCit = GetPatData(model).stateSetBegin();
  int numRates = GetPatData(model).GetNumRates();
  int pveclen = numStates*numRates;
  std::vector<double> categStateProbs(numStates*GetPatData(model).GetNumRates(), 1.0/((double)pveclen));
  /*for (unsigned i = 0; i < numRates; ++i) {
        const double rp = GetPatData(model).GetRateCatProb()[i];
        for (unsigned j = 0; j < numStates; ++j) {
            categStateProbs[i*numStates + j] = rp*(1.0/((double)numStates));
        }
    }*/
  for (; ssCit != GetPatData(model).stateSetEnd(); ssCit++) {
    const int & obsStSet = *ssCit;
    int common = -1;
    while(common > -2){
      ProbForObsStateSet & currNdProbSet = currNdInfo->getForObsStateSet(obsStSet);
      std::vector<double> & currNdProbVec = currNdProbSet.getProbForCommState(common);
      for (int s = 0; s < GetPatData(model).GetNumStates(); s++) {
        //std::cerr << GetPatData(model).GetNumStates() << "\n";
        std::cerr << currNdProbVec[s] << "\n";
        std::cerr << categStateProbs[s] << "\n";
        total += currNdProbVec[s]*categStateProbs[s];
      }
      common = getNextCommStSet(obsStSet, common);
    }
  }
    //std::cerr << "Freeing NodeInfos\n";
    freeNodeInfo(instance, nd, nodeToInfoMap);
    //std::cerr << "Nodes freed\n";
    return total;
    //return currNdInfo;
    //do something else
} // calcUninformativePatterns function end

// returns probability of all patterns being informative
// used for ascertainment bias correction
double totalInformativePatternProb(MTInstance & instance) {
  patClassInitialize(instance);
  //instance.tree.initBrLens();
  double total = 1.0;
  Node * root = instance.tree.GetRoot();

  for (unsigned m = 0; m < instance.numPartitions; m++) {
    unsigned numChars = instance.partMat.GetLeafCharacters(m, 1)->GetCharVec().size();
    // return probability for uniformative pattern for one character
    // should be same for all chars with same number of states
    double uninfCharProb = calcUninformativePatterns(instance, root, 0, m);
    std::cerr << "Partition uninformative prob: " << uninfCharProb << "\n";
    //subtract from total probability to get probability of informative pattern
    double infcharProb = 1.0 - uninfCharProb;
    // raise to power of numChar to get prob of all patterns being informative for that partition
    double partProb = pow(infcharProb, (double) numChars);
    total *= partProb;
    std::cerr << "Partition total: " << partProb << "\n";
  }
  return total;
}



} // namespace
