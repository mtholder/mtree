#include <vector>
#include <cmath>
#include "mt_tree.h"
#include "mt_optimize_branches.h"
#include "mt_tree.h"
#include "mt_data.h"
#include "mt_instance.h"
#include "mt_util.h"

// optimize branch lengths - this is a port of PLL code.
namespace mt {

struct NRSumTable {
    public:
    enum TipCaseEnum tipCase;
    double *x1_start = nullptr;
    double *x2_start = nullptr;
    double *x1_start_asc = nullptr;
    double *x2_start_asc = nullptr;
    unsigned char *tipX1;
    unsigned char *tipX2;
    double *x1_gapColumn = nullptr;
    double *x2_gapColumn = nullptr;
    unsigned int *x1_gap = nullptr;
    unsigned int *x2_gap = nullptr;
    NRSumTable()
        :tipCase(PLL_TIP_TIP),
        x1_start(nullptr),
        x2_start(nullptr),
        x1_start_asc(nullptr),
        x2_start_asc(nullptr),
        tipX1(nullptr),
        tipX2(nullptr),
        x1_gapColumn(nullptr),
        x2_gapColumn(nullptr),
        x1_gap (nullptr),
        x2_gap(nullptr) {
    }
};
// GetNRSumTable = getVects in PLL:
// generic function to get the required pointers to the data associated with the left and right node that define a branch
static NRSumTable GetNRSumTable(MTInstance & instance, std::size_t modelIndex) {
    NRSumTable r;
    const PartitionScoringInfo & psi = instance.GetPartitionScoreInfo(modelIndex);
    const TraversalDescriptor & td = instance.GetTraversalDescriptor();
    const TraversalInfo & ti = td.Top();
        // get the left and right node number of the nodes defining the branch we want to optimize 
    const auto pNumber = td.GetPNumber();
    const auto qNumber = td.GetPNumber();
        // get the index where the ancestral vector is expected to be found
    int p_slot, q_slot;
    const auto mxtips = instance.GetMxTips();
    if (psi.GetUseRecom()) {
        p_slot = ti.slot_p; 
        q_slot = ti.slot_q;
    } else {
        p_slot = pNumber - mxtips - 1;
        q_slot = qNumber - mxtips - 1;
    }
        // switch over the different tip cases again here
    const bool pIsTip = isTip(pNumber, mxtips);
    const bool qIsTip = isTip(qNumber, mxtips);
    const bool doAscBiasCorr = instance.DoAscBiasCorr(modelIndex);

    if (pIsTip || qIsTip) {
        if(!(pIsTip && qIsTip)) {
            r.tipCase = PLL_TIP_INNER;
            const auto yNumber = (qIsTip ? qNumber : pNumber);
            const auto x2Number = (qIsTip ? pNumber : qNumber);
            const auto slotNumber = (qIsTip ? p_slot : q_slot);
            r.tipX1 = psi.GetYArrayPtr(yNumber);
            r.x2_start = psi.GetXArrayPtr(slotNumber);
            if (doAscBiasCorr) {
                r.x2_start_asc = psi.GetAscArrayPtr(x2Number - mxtips - 1); // pointer math... * pr->partitionData[model]->ascOffset
            }
            if (instance.GetSaveMemory()) {
                r.x2_gap = psi.GetGapArrayPtr(x2Number); // pointer math... * pr->partitionData[model]->gapVectorLength
                r.x2_gapColumn = psi.GetGapColumn(x2Number - mxtips - 1); // * states * rateHet];
            }
        } else {
            // note that tip tip should normally not occur since this means that we are trying to optimize 
            // a branch in a two-taxon tree. However, this has been inherited be some RAxML function 
            // that optimized pair-wise distances between all taxa in a tree
            r.tipCase = PLL_TIP_TIP;
            r.tipX1 = psi.GetYArrayPtr(pNumber);
            r.tipX2 = psi.GetYArrayPtr(qNumber);
        }
    } else {
        r.tipCase = PLL_INNER_INNER;
        r.x1_start = psi.GetXArrayPtr(p_slot);
        r.x2_start = psi.GetXArrayPtr(q_slot);
        if (doAscBiasCorr) {
            r.x1_start_asc = psi.GetAscArrayPtr(pNumber - mxtips - 1);
            r.x2_start_asc = psi.GetAscArrayPtr(qNumber - mxtips - 1);
        }
        if (instance.GetSaveMemory()) {
            r.x1_gap = psi.GetGapArrayPtr(pNumber);
            r.x1_gapColumn =psi.GetGapColumn(pNumber - mxtips - 1);
            r.x2_gap = psi.GetGapArrayPtr(qNumber);
            r.x2_gapColumn = psi.GetGapColumn(qNumber - mxtips - 1);
        }
    }
    return r;
}


/** @brief Precompute values (sumtable) from the 2 likelihood vectors of a given branch
 * @warning These precomputations are stored in \a tr->partitionData[model].sumBuffer, which is used by function \a execCore
 * pNumber = tr->td[0].ti[0].pNumber;
 * qNumber = tr->td[0].ti[0].qNumber;
 * @note This function should be called only once at the very beginning of each Newton-Raphson procedure
 * for optimizing barnch lengths. It initially invokes an iterative newview call to get a consistent pair of vectors
 * at the left and the right end of the branch and thereafter invokes the one-time only precomputation of values
 * (sumtable) that can be re-used in each Newton-Raphson iteration. Once this function has been called we can
 * execute the actual NR procedure
 */
void makenewzIterative(MTInstance & instance) {
    int model;
    const auto numberOfPartitions = instance.GetNumPartitions();
    const TraversalDescriptor & td = instance.GetTraversalDescriptor();
    PartitionScoringInfoVec & psiv = instance.GetPartitionScoreInfoVec();
        // call newvieIterative to get the likelihood arrays to the left and right of the branch
    pllNewviewIterative(tr, pr, 1);
        // loop over all partoitions to do the precomputation of the sumTable buffer 
        // This is analogous to the pllNewviewIterative() and pllEvaluateIterative() implementations.
    for(auto model = 0U; model < pr->numberOfPartitions; model++) {
        const auto width = pr->partitionData[model]->width;
        if (td.executeModel[model] && width > 0) {
            const int numStates = psiv.at(model)->numStates;
            NRSumTable sumTable = GetNRSumTable(instance, model);
#           if (!defined(__SSE3) && !defined(__AVX) && !defined(__MIC_NATIVE))
                assert(!tr->saveMemory);
                if(tr->rateHetModel == PLL_CAT) {
                    sumCAT_FLEX(tipCase, width, states, sumTable);
                } else {
                    //sumGAMMA_FLEX_reorder(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                    sumGAMMA_FLEX(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width, states);
                }
#           else
                switch(states) {
                    case 2: /* BINARY */
#                       ifdef __MIC_NATIVE
                            assert(0 && "Binary data model is not implemented on Intel MIC");
#                       else
                            assert(!tr->saveMemory);
                            if (tr->rateHetModel == PLL_CAT)
                                sumCAT_BINARY(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width);
                            else
                                sumGAMMA_BINARY(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width);
#                       endif
                        break;
                    case 4: /* DNA */
#                       ifdef __MIC_NATIVE
                            assert(!tr->saveMemory);
                            assert(tr->rateHetModel == PLL_GAMMA);
                            sumGTRGAMMA_MIC(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width);
#                       else
                            if(tr->rateHetModel == PLL_CAT) {
                                if(tr->saveMemory)
                                  sumCAT_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
                                else
                                  sumCAT(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width);
                                } else {
                                if(tr->saveMemory)
                                  sumGAMMA_GAPPED_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
                                else
                                  sumGAMMA(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width);
                            }
#                       endif
                        break;
                    case 20: /* proteins */
#                       ifdef __MIC_NATIVE
                            assert(!tr->saveMemory);
                            assert(tr->rateHetModel == PLL_GAMMA);
                            if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                                sumGTRGAMMAPROT_LG4_MIC(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector_LG4, tipX1, tipX2, width);
                            else
                                sumGTRGAMMAPROT_MIC(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width);
#                       else
                            if(tr->rateHetModel == PLL_CAT) {
                                if(tr->saveMemory)
                                    sumGTRCATPROT_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
                                else
                                    sumGTRCATPROT(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width);
                            } else {
                                if(tr->saveMemory)
                                    sumGAMMAPROT_GAPPED_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
                                else {
                                    if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                                        sumGAMMAPROT_LG4(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector_LG4, tipX1, tipX2, width);
                                    else
                                        sumGAMMAPROT(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2, width);
                                }
                            }
#                       endif
                        break;
                    default:
                        assert(0);
                }
#           endif
#           if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
                const bool doAscBiasCorr = pr->partitionData[model]->ascBias && tr->threadID == 0;
#           else
                const bool doAscBiasCorr = pr->partitionData[model]->ascBias;
#           endif
            if (doAscBiasCorr) {
                int pNumber = tr->td[0].ti[0].pNumber;
                int qNumber = tr->td[0].ti[0].qNumber;
                const int * const ex1_asc = &pr->partitionData[model]->ascExpVector[(pNumber - tr->mxtips - 1) * states];
                const int * const ex2_asc = &pr->partitionData[model]->ascExpVector[(qNumber - tr->mxtips - 1) * states];
                if (tipCase ==PLL_TIP_INNER) {
                    const double * const exTip = (isTip(pNumber, tr->mxtips) ? ex2_asc : ex1_asc);
                    for (auto i = 0; i < states; i++) {
                        pr->partitionData[model]->ascScaler[i] = pow(PLL_MINLIKELIHOOD, (double) exTip[i]);
                    }
                else if (tipCase == PLL_INNER_INNER) {
                    for (auto i = 0; i < states; i++) {
                        pr->partitionData[model]->ascScaler[i] = pow(PLL_MINLIKELIHOOD, (double) (ex1_asc[i] + ex2_asc[i]));
                    }
                } else {
                  assert(false);
                }
                if (tr->rateHetModel == PLL_CAT)
                    sumCatAsc (tipCase, pr->partitionData[model]->ascSumBuffer, x1_start_asc, x2_start_asc, pr->partitionData[model]->tipVector, states, states);
                else
                    sumGammaAsc(tipCase, pr->partitionData[model]->ascSumBuffer, x1_start_asc, x2_start_asc, pr->partitionData[model]->tipVector, states, states);
            }
        }
    }
}

/* the function below actually implements the iterative Newton-Raphson procedure.
   It is particularly messy and hard to read because for the case of per-partition branch length 
   estimates it needs to keep track of whetehr the Newton Raphson procedure has 
   converged for each partition individually. 
   The rationale for doing it like this is also provided in:
   A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009,

*/
static void topLevelMakenewz(MTInstance & instance,
                             const std::vector<double> & z0, // init branch lengths
                             int _maxiter,
                             std::vector<double> & result) {
    ScoringInfo & si = instance.GetScoringInfo();
    const auto numBranchParts = z0.size();
      // initialize loop convergence variables etc. 
      // maxiter is the maximum number of NR iterations we are going to do before giving up 
    std::vector<double> z = z0;
    std::vector<int> maxiter(numBranchParts, _maxiter);
    std::deque<bool> outerConverged(numBranchParts, false);
    std::deque<bool> curvatOK(numBranchParts, true);
    std::vector<double> coreLogZ(numBranchParts);
    std::vector<double> zstep(numBranchParts);
    std::vector<double> zprev(numBranchParts);
    for (;;) { // check if we are done for partition i or if we need to adapt the branch length again */
        for(auto i = 0U; i < numBranchParts; i++) {
            if (!outerConverged[i]) {
                if (curvatOK[i]) {
                    curvatOK[i] = false;
                    zprev[i] = z[i];
                    zstep[i] = (1.0 - PLL_ZMAX) * z[i] + PLL_ZMIN;
                }
                    // other case, the outer loop hasn't converged but we are trying to approach the maximum from the wrong side 
                enforce_bound(z[i], PLL_ZMIN, PLL_ZMAX);
                coreLogZ[i] = std::log(z[i]);
            }
        SetExecuteMaskByBranchLengthPart(instance, curvatOK, true);
        StoreExecuteMaskInTraversalDescriptor(instance);
        StoreValuesInTraversalDescriptor(instance, coreLogZ);
        std::vector<double> dlnLdlz;
        std::vector<double> d2lnLdlz2;
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
            // if this is the first iteration of NR we will need to first do this one-time call 
            // of maknewzIterative() Note that, only this call requires broadcasting the traversal descriptor,
            // subsequent calls to pllMasterBarrier(PLL_THREAD_MAKENEWZ, tr); will not require this
        if(firstIteration) {
            tr->td[0].traversalHasChanged = true; 
            pllMasterBarrier (tr, pr, PLL_THREAD_MAKENEWZ_FIRST);
            firstIteration = false; 
            tr->td[0].traversalHasChanged = false; 
        } else {
            pllMasterBarrier(tr, pr, PLL_THREAD_MAKENEWZ);
        }
        branchLength_parallelReduce(tr, (double*)dlnLdlz, (double*)d2lnLdlz2, numBranchParts);
#else
            // sequential part, if this is the first newton-raphson implementation,
            // do the precomputations as well, otherwise just execute the computation
            // of the derivatives
        if (firstIteration) {
            makenewzIterative(tr, pr);
            firstIteration = false;
        }
        execCore(tr, pr, dlnLdlz, d2lnLdlz2);
#endif
            // do a NR step, if we are on the correct side of the maximum that's okay, otherwise 
            // shorten branch 
        for(i = 0; i < numBranchParts; i++) {
            if(!outerConverged[i] && !curvatOK[i]) {
                if ((d2lnLdlz2[i] >= 0.0) && (z[i] < PLL_ZMAX)) {
                    zprev[i] = z[i] = 0.37 * z[i] + 0.63; // Bad curvature, shorten branch
                } else {
                    curvatOK[i] = true;
                }
            }
        }
            // do the standard NR step to obrain the next value, depending on the state for eahc partition
        for(i = 0; i < numBranchParts; i++) {
            if (curvatOK[i] && !outerConverged[i]) {
                const double pb = 0.25 * zprev[i] + 0.75;
                if (d2lnLdlz2[i] < 0.0) {
                    double tantmp = -dlnLdlz[i] / d2lnLdlz2[i];
                    if (tantmp < 100) {
                        z[i] *= exp(tantmp);
                    }
                    enforce_min_bound(z[i], PLL_ZMIN);
                    enforce_max_bound(z[i], pb);
                } else {
                    z[i] = pb;
                }
            }
            enforce_max_bound(z[i], PLL_ZMAX);
                // decrement the maximum number of itarations */
            maxiter[i] = maxiter[i] - 1;
                // check if the outer loop has converged */
            if(std::abs(z[i] - zprev[i]) > zstep[i]) {
                   // We should make a more informed decision here, based on the log like improvement
                if (maxiter[i] < -20) {
                    z[i] = z0[i];
                    outerConverged[i] = true;
                } else {
                    outerConverged[i] = false;
                }
            } else {
                outerConverged[i] = true;
            }
        }
        if (all_true(outerConverged)) {
            break;
        }
    }
    SetExecuteMaskForAllPart(instance, true);
    result = z;
}

#if 0
void makenewzGeneric(const unsigned numPartBranchLengths,
                     Arc edge,
                     int maxiter,
                     double *result,
                     bool checkPartition) {
    bool p_recom = false; /* if one of was missing, we will need to force recomputation */
    bool q_recom = false;
    /* the first entry of the traversal descriptor stores the node pair that defines the branch */
    tr->td[0].ti[0].pNumber = p->number;
    tr->td[0].ti[0].qNumber = q->number;
    for(auto i = 0U; i < numPartBranchLengths; i++) {
        //originalExecute[i] = pr->partitionData[i]->executeModel;
        tr->td[0].ti[0].qz[i] = z0[i];
        if(mask) {
            pr->partitionData[i]->executeModel = !tr->partitionConverged[i]
        }
    }
    if (tr->useRecom)
    {
    int
    slot = -1;
    //count = 0;

    /* Ensure p and q get a unpinnable slot in physical memory */
    if(!isTip(q->number, tr->mxtips))
    {
    q_recom = getxVector(tr->rvec, q->number, &slot, tr->mxtips);
    tr->td[0].ti[0].slot_q = slot;
    }
    if(!isTip(p->number, tr->mxtips))
    {
    p_recom = getxVector(tr->rvec, p->number, &slot, tr->mxtips);
    tr->td[0].ti[0].slot_p = slot;
    }
    }


    /* compute the traversal descriptor of the likelihood vectors that need to be re-computed
    first in makenewzIterative */

    tr->td[0].count = 1;

    if(p_recom || needsRecomp(tr->useRecom, tr->rvec, p, tr->mxtips))
    computeTraversal(tr, p, PLL_TRUE, numBranches);

    if(q_recom || needsRecomp(tr->useRecom, tr->rvec, q, tr->mxtips))
    computeTraversal(tr, q, PLL_TRUE, numBranches);

    /* call the Newton-Raphson procedure */

    topLevelMakenewz(tr, pr, z0, maxiter, result);

    /* Mark node as unpinnable */
    if(tr->useRecom)
    {
    unpinNode(tr->rvec, p->number, tr->mxtips);
    unpinNode(tr->rvec, q->number, tr->mxtips);
    }

    /* fix eceuteModel this seems to be a bit redundant with topLevelMakenewz */

    for(i = 0; i < numBranches; i++)
    pr->partitionData[i]->executeModel = PLL_TRUE;
}
#endif

double optimizeAllLengthsForOneEdge(Arc edge, MTInstance & mtInstance) {
    const unsigned numPartBranchLengths = mtInstance.GetNumPartitions();
    const double brLenDiffThreshold = mtInstance.optSettings.brLenDiffThreshold;
    auto elvp = edge.GetEdgeLenVec();
    assert(elvp != nullptr);
    assert(elvp->size() >= numPartBranchLengths);
    std::vector<double> edgeLenCopy{*elvp};
    int newzpercycle = 1;
    makenewzGeneric(numPartBranchLengths, edge, newzpercycle, &edgeLenCopy[0]);
    std::vector<bool> smoothed = mtInstance.optSettings.partitionSmoothed;
    for(auto i = 0U; i < numPartBranchLengths; i++) {
        if (mtInstance.optSettings.partitionConverged[i]) {
            if(abs(edgeLenCopy[i] - (*elvp)[i]) > brLenDiffThreshold) {
                smoothed[i] = 0;
            }
            (*elvp)[i] = edgeLenCopy[i];
        }
    }
    mtInstance.optSettings.partitionSmoothed = smoothed;
}

} // namespace mt

#if 0
// optimize branch lengths - this is a port of PLL code.
namespace mt {

/// smoothTree in PLL

static bool allSmoothed (MTInstance &instance, int numBranches) {
  bool result = true;

  for(int i = 0; i < numBranches; i++)
  {
    if(instance.optSettings.partitionSmoothed[i] == false)
      result = false;
    else
      instance.optSettings.partitionConverged[i] = true;
  }

  return result;
}
}

double optimizeAllBranchLengthsForAllPartitions(MTInstance &instance) {

    const unsigned numPartBranchLengths = instance.partMat.GetNumPartitions();
    instance.optSettings.partitionConverged.assign(numPartBranchLengths, 0);
    unsigned maxLoops = instance.optSettings.maxIterBrLenSmoothing;
    while (--maxLoops >= 0) {
        instance.optSettings.partitionSmoothed.assign(numPartBranchLengths, 1);
        sweepOverTreeOptimizeAllBranchLengthsForAllPartitions(instance);
        if (allSmoothed(instance, numPartBranchLengths))
            break;
    }
    instance.settings.partitionConverged.assign(numPartBranchLengths, 0);
}


template<typename T>
class BeforeAfterIter {
    BeforeAfterIter(Node *, (*before)(Node *, T), (*after)(Node *, T), T blob);
};

class BeforeAfterIter {
    BeforeAfterIter(Node *, (*before)(Node *, void *), (*after)(Node *, void *), void * blob);
};

/** sweepOverTreeOptimizeAllBranchLengthsForAllPartitions = smooth in PLL
*/
void sweepOverTreeOptimizeAllBranchLengthsForAllPartitions (MTInstance &instance)
{
    const unsigned numPartBranchLengths = instance.partMat.GetNumPartitions();
    Tree & tree = instance.tree;
    Node * root = tree.GetRoot();
    void * blob = &(instance);
    BeforeAfterIter edgeIter(root, optimizeAllLengthsForOneEdgeHook, updatCLAForExit, blob);
    while (edgeIter.get() != nullptr) {
        edgeIter.advance();
    }
}

void updatePartials(MTInstance &instance, Node * node) {

}

bool updatCLAForExit(Node * node, void * mtInstance) {
    if(numBranches > 1 && !tr->useRecom)
      updatePartials(mtInstance, node);
    else
      updatePartials(mtInstance, node);
}

//pll optimizeAllLengthsForOneEdge update
bool optimizeAllLengthsForOneEdgeHook(Node * node, void * mtInstance) {
    MTInstance * instance = (MTInstance *)mtInstance;
    Arc arc = (node, node->parent);
    return optimizeAllLengthsForOneEdge(arc, *instance);
}


/*
void update(MTInstance &instance, Arc edge)
{
  double z[PLL_NUM_BRANCHES], z0[PLL_NUM_BRANCHES];
  const unsigned numBranches = instance.partMat.GetNumPartitions();


  q = p->back;

  for(i = 0; i < numBranches; i++)
    z0[i] = q->z[i];

  if(numBranches > 1)
    makenewzGeneric(tr, pr, p, q, z0, PLL_NEWZPERCYCLE, z, PLL_TRUE);
  else
    makenewzGeneric(tr, pr, p, q, z0, PLL_NEWZPERCYCLE, z, PLL_FALSE);

  for(i = 0; i < numBranches; i++)
  {
    if(!tr->partitionConverged[i])
    {
      if(PLL_ABS(z[i] - z0[i]) > PLL_DELTAZ)
      {
        tr->partitionSmoothed[i] = PLL_FALSE;
      }

      p->z[i] = q->z[i] = z[i];
    }
  }

}
*/


#endif

#if 0
// makenewzGeneric in PLL (makenewzgenericspecial.c)
/** @brief Optimize branch length value(s) of a given branch with the Newton-Raphtson procedure
 *
 * @warning A given branch may have one or several branch length values (up to PLL_NUM_BRANCHES), usually the later refers to partition-specific branch length values. Thus z0 and result represent collections rather than double values. The number of branch length values is given by \a tr->numBranches
 *
 * @param tr
 *   Library instance
 *
 * @param p
 *   One node that defines the branch (p->z)
 *
 * @param q
 *   The other node side of the branch (usually p->back), but the branch length can be estimated even if p and q are
 *   not connected, e.g. before the insertion of a subtree.
 *
 * @param z0
 *   Initial branch length value(s) for the given branch \a p->z
 *
 * @param maxiter
 *   Maximum number of iterations in the Newton-Raphson procedure
 *
 * @param result
 *   Resulting branch length value(s) for the given branch \a p->z
 *
 * @param mask
 *   Specifies if a mask to track partition convergence (\a tr->partitionConverged) is being used.
 *
 * @sa typical values for \a maxiter are constants \a iterations and \a PLL_NEWZPERCYCLE
 * @note Requirement: q->z == p->z
 */
void makenewzGeneric(pllInstance *tr, partitionList * pr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, boolean mask)
{
  int i;
  //boolean originalExecute[PLL_NUM_BRANCHES];
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  boolean
    p_recom = PLL_FALSE, /* if one of was missing, we will need to force recomputation */
    q_recom = PLL_FALSE;

  /* the first entry of the traversal descriptor stores the node pair that defines
     the branch */

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(i = 0; i < numBranches; i++)
  {
    //originalExecute[i] = pr->partitionData[i]->executeModel;
    tr->td[0].ti[0].qz[i] = z0[i];
    if(mask)
    {
      if (tr->partitionConverged[i])
        pr->partitionData[i]->executeModel = PLL_FALSE;
      else
        pr->partitionData[i]->executeModel = PLL_TRUE;
    }
  }
  if (tr->useRecom)
  {
    int
      slot = -1;
      //count = 0;

    /* Ensure p and q get a unpinnable slot in physical memory */
    if(!isTip(q->number, tr->mxtips))
    {
      q_recom = getxVector(tr->rvec, q->number, &slot, tr->mxtips);
      tr->td[0].ti[0].slot_q = slot;
    }
    if(!isTip(p->number, tr->mxtips))
    {
      p_recom = getxVector(tr->rvec, p->number, &slot, tr->mxtips);
      tr->td[0].ti[0].slot_p = slot;
    }
  }


  /* compute the traversal descriptor of the likelihood vectors that need to be re-computed
     first in makenewzIterative */

  tr->td[0].count = 1;

  if(p_recom || needsRecomp(tr->useRecom, tr->rvec, p, tr->mxtips))
    computeTraversal(tr, p, PLL_TRUE, numBranches);

  if(q_recom || needsRecomp(tr->useRecom, tr->rvec, q, tr->mxtips))
    computeTraversal(tr, q, PLL_TRUE, numBranches);

  /* call the Newton-Raphson procedure */

  topLevelMakenewz(tr, pr, z0, maxiter, result);

  /* Mark node as unpinnable */
  if(tr->useRecom)
  {
    unpinNode(tr->rvec, p->number, tr->mxtips);
    unpinNode(tr->rvec, q->number, tr->mxtips);
  }

  /* fix eceuteModel this seems to be a bit redundant with topLevelMakenewz */

  for(i = 0; i < numBranches; i++)
    pr->partitionData[i]->executeModel = PLL_TRUE;
}

/* the function below actually implements the iterative Newton-Raphson procedure.
   It is particularly messy and hard to read because for the case of per-partition branch length
   estimates it needs to keep track of whether the Newton Raphson procedure has
   converged for each partition individually.

   The rational efor doing it like this is also provided in:


   A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009,

*/


void topLevelMakenewz(MTInstance *tr, partitionList * pr, double *z0, int _maxiter, double *result)
{
  double z[PLL_NUM_BRANCHES], zprev[PLL_NUM_BRANCHES], zstep[PLL_NUM_BRANCHES];
  volatile double dlnLdlz[PLL_NUM_BRANCHES], d2lnLdlz2[PLL_NUM_BRANCHES];
  int i, maxiter[PLL_NUM_BRANCHES], model;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;
  boolean firstIteration = PLL_TRUE;
  boolean outerConverged[PLL_NUM_BRANCHES];
  boolean loopConverged;

  /* figure out if this is on a per partition basis or jointly across all partitions */
  /* initialize loop convergence variables etc.
     maxiter is the maximum number of NR iterations we are going to do before giving up */

  for(i = 0; i < numBranches; i++)
  {
    z[i] = z0[i];
    maxiter[i] = _maxiter;
    outerConverged[i] = PLL_FALSE;
    tr->curvatOK[i] = PLL_TRUE;
  }


  /* nested do while loops of Newton-Raphson */
  do
  {
    /* check if we ar done for partition i or if we need to adapt the branch length again */
    for(i = 0; i < numBranches; i++)
    {
      if(outerConverged[i] == PLL_FALSE && tr->curvatOK[i] == PLL_TRUE)
      {
        tr->curvatOK[i] = FALSE;

        zprev[i] = z[i];

        zstep[i] = (1.0 - PLL_ZMAX) * z[i] + PLL_ZMIN;
      }
    }

    for(i = 0; i < numBranches; i++)
    {
      /* other case, the outer loop hasn't converged but we are trying to approach
         the maximum from the wrong side */

      if(outerConverged[i] == PLL_FALSE && tr->curvatOK[i] == PLL_FALSE)
      {
        double lz;

        if (z[i] < PLL_ZMIN) z[i] = PLL_ZMIN;
        else if (z[i] > PLL_ZMAX) z[i] = PLL_ZMAX;
        lz = log(z[i]);

        tr->coreLZ[i] = lz;
      }
    }


    /* set the execution mask */

    if(numBranches > 1)
    {
      for(model = 0; model < pr->numberOfPartitions; model++)
      {
        if(pr->partitionData[model]->executeModel)
          pr->partitionData[model]->executeModel = !tr->curvatOK[model];

      }
    }
    else
    {
      for(model = 0; model < pr->numberOfPartitions; model++)
        pr->partitionData[model]->executeModel = !tr->curvatOK[0];
    }


    /* store it in traversal descriptor */

    storeExecuteMaskInTraversalDescriptor(tr, pr);

    /* store the new branch length values to be tested in traversal descriptor */

    storeValuesInTraversalDescriptor(tr, pr, &(tr->coreLZ[0]));

//#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

    /* if this is the first iteration of NR we will need to first do this one-time call
       of maknewzIterative() Note that, only this call requires broadcasting the traversal descriptor,
       subsequent calls to pllMasterBarrier(PLL_THREAD_MAKENEWZ, tr); will not require this
       */

    if(firstIteration)
      {
        tr->td[0].traversalHasChanged = PLL_TRUE;
        pllMasterBarrier (tr, pr, PLL_THREAD_MAKENEWZ_FIRST);
        firstIteration = PLL_FALSE;
        tr->td[0].traversalHasChanged = PLL_FALSE;
      }
    else
      pllMasterBarrier(tr, pr, PLL_THREAD_MAKENEWZ);
    branchLength_parallelReduce(tr, (double*)dlnLdlz, (double*)d2lnLdlz2, numBranches);
    printf(" makenewz after branchLength_parallelReduce dlnLdlz = %lf d2lnLdlz2 = %lf\n", dlnLdlz, d2lnLdlz2);
//#else
    /* sequential part, if this is the first newton-raphson implementation,
       do the precomputations as well, otherwise just execute the computation
       of the derivatives */
    if(firstIteration)
      {
        makenewzIterative(tr, pr);
        firstIteration = FALSE;
      }
    execCore(tr, pr, dlnLdlz, d2lnLdlz2);
    printf(" makenewz after execCore dlnLdlz = %lf d2lnLdlz2 = %lf\n", dlnLdlz, d2lnLdlz2);
#endif

#if 0
    /* do a NR step, if we are on the correct side of the maximum that's okay, otherwise
       shorten branch */

    for(i = 0; i < numBranches; i++)
    {
      if(outerConverged[i] == PLL_FALSE && tr->curvatOK[i] == PLL_FALSE)
      {
        if ((d2lnLdlz2[i] >= 0.0) && (z[i] < PLL_ZMAX))
          zprev[i] = z[i] = 0.37 * z[i] + 0.63; /* Bad curvature, shorten branch */
        else
          tr->curvatOK[i] = PLL_TRUE;
      }
    }

    /* do the standard NR step to obrain the next value, depending on the state for eahc partition */

    for(i = 0; i < numBranches; i++)
    {
      if(tr->curvatOK[i] == PLL_TRUE && outerConverged[i] == PLL_FALSE)
      {
        if (d2lnLdlz2[i] < 0.0)
        {
          double tantmp = -dlnLdlz[i] / d2lnLdlz2[i];
          if (tantmp < 100)
          {
            z[i] *= exp(tantmp);
            if (z[i] < PLL_ZMIN)
              z[i] = PLL_ZMIN;

            if (z[i] > 0.25 * zprev[i] + 0.75)
              z[i] = 0.25 * zprev[i] + 0.75;
          }
          else
            z[i] = 0.25 * zprev[i] + 0.75;
        }
        if (z[i] > PLL_ZMAX)
            z[i] = PLL_ZMAX;

        /* decrement the maximum number of itarations */

        maxiter[i] = maxiter[i] - 1;

        /* check if the outer loop has converged */

        //old code below commented out, integrated new PRELIMINARY BUG FIX !
        //this needs further work at some point!

        /*
        if(maxiter[i] > 0 && (PLL_ABS(z[i] - zprev[i]) > zstep[i]))
          outerConverged[i] = PLL_FALSE;
        else
          outerConverged[i] = PLL_TRUE;
        */

        if((PLL_ABS(z[i] - zprev[i]) > zstep[i]))
         {
           /* We should make a more informed decision here,
              based on the log like improvement */

           if(maxiter[i] < -20)
            {
              z[i] = z0[i];
              outerConverged[i] = PLL_TRUE;
            }
           else
             outerConverged[i] = PLL_FALSE;
         }
        else
          outerConverged[i] = PLL_TRUE;
      }
    }

    /* check if the loop has converged for all partitions */

    loopConverged = PLL_TRUE;
    for(i = 0; i < numBranches; i++)
      loopConverged = loopConverged && outerConverged[i];
  }
  while (!loopConverged);


  /* reset partition execution mask */

  for(model = 0; model < pr->numberOfPartitions; model++)
    pr->partitionData[model]->executeModel = PLL_TRUE;

  /* copy the new branches in the result array of branches.
     if we don't do a per partition estimate of
     branches this will only set result[0]
     */

  for(i = 0; i < numBranches; i++)
    result[i] = z[i];
}



}// namespace
#endif
