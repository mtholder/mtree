#include <vector>
#if 0
#include "mt_tree.h"
#include "mt_optimize_branches.h"
#include "mt_instance.h"
#include "mt_tree.h"
#include "mt_data.h"


// optimize branch lengths - this is a port of PLL code.
namespace mt {

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

void optimizeAllBranchLengthsForAllPartitions(mt::MTInstance &instance) {

    const unsigned numPartBranchLengths = instance.partMat.GetNumPartitions();
    instance.optSettings.partitionConverged.assign(numPartBranchLengths, false);
    unsigned maxLoops = instance.optSettings.maxIterBrLenSmoothing;
    while (--maxLoops >= 0) {
        instance.optSettings.partitionSmoothed.assign(numPartBranchLengths, 1);
        sweepOverTreeOptimizeAllBranchLengthsForAllPartitions(instance);
        if (allSmoothed(instance, numPartBranchLengths))
            break;
    }
    instance.settings.partitionConverged.assign(numPartBranchLengths, 0);
}

/*
template<typename T>
class BeforeAfterIter {
    BeforeAfterIter(Node *, (*before)(Node *, T), (*after)(Node *, T), T blob);
};
*/

typedef Node * (*IterFunc)(Node *, void *);

class BeforeAfterIter {
    public:
      BeforeAfterIter(Node * start, IterFunc bef, IterFunc aft, void * blob) {
        startNd(start),
        instptr(blob),
        before(bef),
        after(aft)
      };
      Node * get() {
        return this->start;
      };
      void advance() {
        start = start->leftChild;
      }
    private:
      Node * start;
      void * instptr;
      IterFunc before;
      IterFunc after;
};


/**  sweepOverTreeOptimizeAllBranchLengthsForAllPartitions = smooth in PLL
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

bool optimizeAllLengthsForOneEdge(Arc edge, MTInstance & mtInstance) {
    const unsigned numPartBranchLengths = instance.GetNumPartitions();
    const double brLenDiffThreshold = mtInstance.optSettings.brLenDiffThresh;
    vector<double> origEdgeLen(edge.GetEdgeLengthValue());
    vector<double> edgeLenCopy(origEdgeLen);
    int newzpercycle = 1;
    makenewzGeneric(mtInstance, edge, newzpercycle, &edgeLenCopy[0]);
    vector<char> smoothed = instance.optSettings.partitionSmoothed;
    for(auto i = 0U; i < numPartBranchLengths; i++) {
        if (instance.settings.partitionConverged[i]) {
            if(abs(edgeLenCopy[i] - origEdgeLen[i]) > brLenDiffThresh) {
                smoothed[i] = 0;
            }
            edgeLengthPtr[i] = edgeLenCopy[i];
        }
    }
    instance.optSettings.partitionSmoothed = smoothed;
  return TRUE;
}


void update(MTInstance &instance, Arc edge)
{
  double   z, z0;
//  const unsigned numBranches = instance.partMat.GetNumPartitions();

  p = edge.toNode;
  q = edge.fromNode;

//  for(i = 0; i < numBranches; i++)
    z0 = q.GetEdgeLen();

  if(numBranches > 1)
    makenewzGeneric());
  else
    makenewzGeneric();

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

// makenewzIterative in PLL
void nrPrecomp(MTInstance *instance, Arc edge)
{
  int model, tipCase;

  double
    *x1_start     = NULL,
    *x2_start     = NULL,
    *x1_start_asc = NULL,
    *x2_start_asc = NULL;


  unsigned char
    *tipX1,
    *tipX2;

  double
    *x1_gapColumn = (double*)NULL,
    *x2_gapColumn = (double*)NULL;

  unsigned int
    *x1_gap = (unsigned int*)NULL,
    *x2_gap = (unsigned int*)NULL;

  /* call newvieIterative to get the likelihood arrays to the left and right of the branch */

  pllNewviewIterative(tr, pr, 1);


  /*
     loop over all partoitions to do the precomputation of the sumTable buffer
     This is analogous to the pllNewviewIterative() and pllEvaluateIterative()
     implementations.
     */

  for(model = 0; model < pr->numberOfPartitions; model++)
  {
    int
      width = pr->partitionData[model]->width;

    if(tr->td[0].executeModel[model] && width > 0)
    {
      int
        states = pr->partitionData[model]->states;


      getVects(tr, pr, &tipX1, &tipX2, &x1_start, &x2_start, &tipCase, model, &x1_gapColumn, &x2_gapColumn, &x1_gap, &x2_gap, &x1_start_asc, &x2_start_asc);

#if (!defined(__SSE3) && !defined(__AVX) && !defined(__MIC_NATIVE))
      assert(!tr->saveMemory);
      if(tr->rateHetModel == PLL_CAT)
        sumCAT_FLEX(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
            width, states);
      else
        //sumGAMMA_FLEX_reorder(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
          sumGAMMA_FLEX(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
            width, states);
#else
      switch(states)
      {
      case 2: /* BINARY */
#ifdef __MIC_NATIVE
          assert(0 && "Binary data model is not implemented on Intel MIC");
#else
          assert(!tr->saveMemory);
          if (tr->rateHetModel == PLL_CAT)
            sumCAT_BINARY(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                          width);
          else
            sumGAMMA_BINARY(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                            width);
#endif
          break;
      case 4: /* DNA */
#ifdef __MIC_NATIVE
      assert(!tr->saveMemory);
      assert(tr->rateHetModel == PLL_GAMMA);

      sumGTRGAMMA_MIC(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
          width);
#else
          if(tr->rateHetModel == PLL_CAT)
          {
            if(tr->saveMemory)
              sumCAT_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else
              sumCAT(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width);
          }
          else
          {
            if(tr->saveMemory)
              sumGAMMA_GAPPED_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else
              sumGAMMA(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width);
          }
#endif
          break;
        case 20: /* proteins */
#ifdef __MIC_NATIVE
          assert(!tr->saveMemory);
          assert(tr->rateHetModel == PLL_GAMMA);

              if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                          sumGTRGAMMAPROT_LG4_MIC(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector_LG4, tipX1, tipX2,
                                  width);
              else
                          sumGTRGAMMAPROT_MIC(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                                  width);
#else

            if(tr->rateHetModel == PLL_CAT)
          {
            if(tr->saveMemory)
              sumGTRCATPROT_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector,
                  tipX1, tipX2, width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else
              sumGTRCATPROT(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector,
                  tipX1, tipX2, width);
          }
          else
          {

            if(tr->saveMemory)
              sumGAMMAPROT_GAPPED_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
              else
                    {
                      if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                        sumGAMMAPROT_LG4(tipCase,  pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector_LG4,
                                         tipX1, tipX2, width);
            else
              sumGAMMAPROT(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector,
                  tipX1, tipX2, width);
                    }
          }
#endif
          break;
        default:
          assert(0);
      }
#endif
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      if (pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
      if (pr->partitionData[model]->ascBias)
#endif
       {
            int pNumber = tr->td[0].ti[0].pNumber, qNumber =
                    tr->td[0].ti[0].qNumber, i, *ex1_asc =
                    &pr->partitionData[model]->ascExpVector[(pNumber
                            - tr->mxtips - 1) * states], *ex2_asc =
                    &pr->partitionData[model]->ascExpVector[(qNumber
                            - tr->mxtips - 1) * states];
            switch (tipCase)
            {
            case PLL_TIP_TIP:
                assert(0);
                break;
            case PLL_TIP_INNER:
                if (isTip(pNumber, tr->mxtips))
                {
                    for (i = 0; i < states; i++)
                        pr->partitionData[model]->ascScaler[i] = pow(
                                PLL_MINLIKELIHOOD, (double) ex2_asc[i]);
                }
                else
                {
                    for (i = 0; i < states; i++)
                        pr->partitionData[model]->ascScaler[i] = pow(
                                PLL_MINLIKELIHOOD, (double) ex1_asc[i]);
                }
                break;
            case PLL_INNER_INNER:
                for (i = 0; i < states; i++)
                    pr->partitionData[model]->ascScaler[i] = pow(
                            PLL_MINLIKELIHOOD,
                            (double) (ex1_asc[i] + ex2_asc[i]));
                break;
            default:
                assert(0);
            }
         if (tr->rateHetModel == PLL_CAT)
           sumCatAsc  (tipCase, pr->partitionData[model]->ascSumBuffer, x1_start_asc, x2_start_asc, pr->partitionData[model]->tipVector, states, states);
         else
           sumGammaAsc(tipCase, pr->partitionData[model]->ascSumBuffer, x1_start_asc, x2_start_asc, pr->partitionData[model]->tipVector, states, states);
       }
    }
  }
}

// based on makenewzGeneric in PLL
void optBrLen(MTInstance &instance, Arc edge, int maxiter, double *result)
{
  int i;
  //boolean originalExecute[PLL_NUM_BRANCHES];
  const unsigned numPartBranchLengths = instance.partMat.GetNumPartitions();

  bool
    p_recom = false, /* if one of was missing, we will need to force recomputation */
    q_recom = false;

  /* the first entry of the traversal descriptor stores the node pair that defines
     the branch */


  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(i = 0; i < numBranches; i++)
  {
    //originalExecute[i] =  pr->partitionData[i]->executeModel;
    tr->td[0].ti[0].qz[i] =  z0[i];
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


  /* compute the traversal descriptor of the likelihood vectors  that need to be re-computed
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
#if 0
/// topLevelMakenewz in PLL
void newtRaph (MTInstance *instance, Arc edge, double *z0, int maxiter, double *result)
{
  //double   z[50], zprev[50], zstep[50];
  //volatile double  dlnLdlz[50], d2lnLdlz2[50];

  double z, zprev, zstep;
  volatile double dlnLdlz, l2lnLdlz2;

  int i, maxiter[50], model;

  //int numBranches = instance.partMat.GetNumPartitions();

  bool firstIteration = TRUE;

  //bool outerConverged[PLL_NUM_BRANCHES];
  bool outerConverged = FALSE;

  bool loopConverged = FALSE;

  /* nested do while loops of Newton-Raphson */
  do
  {
    /* check if we are done or if we need to adapt the branch length again */
//    if(!outerConverged && instance.curvatOK)
//    {
//        instance.curvatOK = FALSE;

//        zprev[i] = z[i];

//        zstep[i] = (1.0 - PLL_ZMAX) * z[i] + PLL_ZMIN;
//    }

      /* other case, the outer loop hasn't converged but we are trying to approach
         the maximum from the wrong side */

//      if(outerConverged[i] == PLL_FALSE && tr->curvatOK[i] == PLL_FALSE)
//      {
//        double lz;

//        if (z[i] < PLL_ZMIN) z[i] = PLL_ZMIN;
//        else if (z[i] > PLL_ZMAX) z[i] = PLL_ZMAX;
//        lz    = log(z[i]);

//        tr->coreLZ[i] = lz;
//      }

    /* set the execution mask */

//    if(numBranches > 1)
//    {
//      for(model = 0; model < pr->numberOfPartitions; model++)
//      {
//        if(pr->partitionData[model]->executeModel)
//          pr->partitionData[model]->executeModel = !tr->curvatOK[model];

//      }
//    }
//    else
//    {
//      for(model = 0; model < pr->numberOfPartitions; model++)
//        pr->partitionData[model]->executeModel = !tr->curvatOK[0];
//    }


    /* store it in traversal descriptor */

//    storeExecuteMaskInTraversalDescriptor(tr, pr);

    /* store the new branch length values to be tested in traversal descriptor */

//    storeValuesInTraversalDescriptor(tr, pr, &(tr->coreLZ[0]));

    /* if this is the first iteration of NR we will need to first do this one-time call
       of maknewzIterative() Note that, only this call requires broadcasting the traversal descriptor,
       subsequent calls to pllMasterBarrier(PLL_THREAD_MAKENEWZ, tr); will not require this
       */

//    if(firstIteration)
//      {
//        tr->td[0].traversalHasChanged = PLL_TRUE;
//        pllMasterBarrier (tr, pr, PLL_THREAD_MAKENEWZ_FIRST);
//        firstIteration = PLL_FALSE;
//        tr->td[0].traversalHasChanged = PLL_FALSE;
//      }
//    else
//      pllMasterBarrier(tr, pr, PLL_THREAD_MAKENEWZ);
//    branchLength_parallelReduce(tr, (double*)dlnLdlz, (double*)d2lnLdlz2, numBranches);
//    printf(" makenewz after branchLength_parallelReduce dlnLdlz = %lf d2lnLdlz2 = %lf\n", dlnLdlz, d2lnLdlz2);

    /* sequential part, if this is the first newton-raphson implementation,
       do the precomputations as well, otherwise just execute the computation
       of the derivatives */
    if(firstIteration)
      {
        nrPrecomp(instance, edge);
        firstIteration = FALSE;
      }
    execCore(tr, pr, dlnLdlz, d2lnLdlz2);
    printf(" makenewz after execCore dlnLdlz = %lf d2lnLdlz2 = %lf\n", dlnLdlz, d2lnLdlz2);

    /* do a NR step, if we are on the correct side of the maximum that's okay, otherwise
       shorten branch */

    for(i = 0; i < numBranches; i++)
    {
      if(outerConverged[i] == PLL_FALSE && tr->curvatOK[i] == PLL_FALSE)
      {
        if ((d2lnLdlz2[i] >= 0.0) && (z[i] < PLL_ZMAX))
          zprev[i] = z[i] = 0.37 * z[i] + 0.63;  /*  Bad curvature, shorten branch */
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
    loopConverged = loopConverged && outerConverged[i];
  }
  while (!loopConverged);


  /* reset  partition execution mask */

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

