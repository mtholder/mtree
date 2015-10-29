#include "mt_instance.h"
namespace mt {

MTInstance::MTInstance(const std::map<unsigned, std::vector< std::vector<mt::char_state_t> > > & ns2RawMat,
                       const std::map<unsigned, mt::CharStateToPrimitiveInd > & ns2cs2pi,
                       const std::vector<double> &patternWts,
                       const std::vector<std::size_t> &orig2compressed,
                       const std::map<unsigned, std::set<unsigned> > & numStates2PatternIndexSet,
                       const NxsSimpleTree & nxsTree, // should be adapted in main
                       const ModelDescription & md,
                       const std::map<unsigned, std::size_t> & numStates2NumBogusChar)
    :partMat(static_cast<unsigned>(ns2RawMat.begin()->second.size()),
             patternWts,
             orig2compressed,
             numStates2PatternIndexSet,
             numStates2NumBogusChar),
    tree(2*static_cast<unsigned>(ns2RawMat.begin()->second.size()) - 1, static_cast<unsigned>(ns2RawMat.begin()->second.size())),
    HasSearchConverged(false),
    curvatOK(true),
    numPartitions(static_cast<unsigned>(numStates2PatternIndexSet.size())),
    likelihoods(numStates2PatternIndexSet.size(),0.0),
    dirtyFlags(numStates2PatternIndexSet.size(), true),
    models(numStates2PatternIndexSet.size(), nullptr) {
    const unsigned numTaxa = static_cast<unsigned>(ns2RawMat.begin()->second.size());
    const unsigned numRateCats = 4; // TEMP: change to customizable value in INI

    std::map<const NxsSimpleNode *, unsigned> ncl2nodeNumber;
    std::vector<const NxsSimpleNode *> pre = nxsTree.GetPreorderTraversal();
    unsigned internalIndex = numTaxa;
    for (std::vector<const NxsSimpleNode *>::iterator ndIt = pre.begin(); ndIt != pre.end(); ++ndIt) {
        const NxsSimpleNode *nd = *ndIt;
        unsigned num;
        if (nd->IsTip()) {
            num = nd->GetTaxonIndex();
        } else {
            num = internalIndex++;
        }
        mt::Node * treeNode = tree.GetNode(num);
        const NxsSimpleNode * par = nd->GetEdgeToParent().GetParent();
        if (par == nullptr) {
            tree.SetRoot(treeNode);
        } else {
            assert (ncl2nodeNumber.find(par) != ncl2nodeNumber.end());
            unsigned parNodeNumber = ncl2nodeNumber[par];
            mt::Node * parNode = tree.GetNode(parNodeNumber);
            parNode->AddChild(treeNode, nd->GetEdgeToParent().GetDblEdgeLen());
        }
        ncl2nodeNumber[nd] = num;
    }

    unsigned partIndex = 0;
    const auto abm = md.GetAscBiasMode();
    unsigned numNodes = 2 * numTaxa - 1;
    //BitFieldMatrix bMat;
    for (auto ns2pi: numStates2PatternIndexSet) {
        const unsigned currPartNumStates = ns2pi.first;
        const auto & rawMat = ns2RawMat.at(currPartNumStates);
        mt::CharModel *m = nullptr;
        if (abm == mt::ModelDescription::NO_ASC_BIAS) {
            m = new mt::MkCharModel(currPartNumStates, numRateCats);
        } else if (abm == mt::ModelDescription::VAR_ONLY_NO_MISSING_ASC_BIAS) {
            m = new mt::MkVarNoMissingAscCharModel(currPartNumStates, numRateCats);
        } else if (abm == mt::ModelDescription::VAR_ONLY_MISSING_ASC_BIAS) {
            m = new mt::MkVarMissingAscCharModel(currPartNumStates, numRateCats);
        } else if (abm == mt::ModelDescription::PARS_ONLY_NO_MISSING_ASC_BIAS) {
            m = new mt::MkParsInfNoMissingModel(currPartNumStates, numRateCats);
        } else if (abm == mt::ModelDescription::PARS_ONLY_MISSING_ASC_BIAS) {
            m = new mt::MkParsInfMissingModel(currPartNumStates, numRateCats);
        } else {
            assert(false);
            std::cerr << "Unrecognoized ASC BIAS MODE\n";
            std::exit(1);
        }
        //m->alphabet = convertToBitFieldMatrix(*charsBlock, bMat);
        models[partIndex] = m;
        const auto & cs2pi = ns2cs2pi.at(currPartNumStates);
        partMat.fillPartition(partIndex, rawMat, &cs2pi);
        for (auto li = 0U; li < numTaxa; ++li) {
            mt::Node * leaf = tree.GetLeaf(li);
            assert(leaf);
            //std::cerr << "Number of chars in this partition: " << rawMat[1].size() << '\n';
            leaf->SetData(partIndex, (void *) partMat.GetLeafCharacters(partIndex, li));
            leaf->SetWork(partIndex, (void *) new mt::LeafWork(rawMat[1].size(),
                                                               static_cast<unsigned>(cs2pi.size()),
                                                               currPartNumStates,
                                                               numRateCats));
        }

        for (auto li = numTaxa; li < numNodes; ++ li) {
            mt::Node * nd = tree.GetNode(li);
            assert(nd);
            nd->SetWork(partIndex, (void *) new mt::InternalNodeWork(rawMat[1].size(),
                                                                     currPartNumStates,
                                                                     numRateCats));
        }

        ++partIndex;
    }
    //exit(1);
}

MTInstance::~MTInstance() {
    for (auto m : models) {
        delete m;
    }
    models.clear();
}

} // namespace
