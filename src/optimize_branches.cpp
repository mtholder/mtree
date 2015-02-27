#include "mt_optimize_branches.h"
// optimize branch lengths - this is a port of PLL code.
namespace mt {

/// smoothTree in PLL

void optimizeAllBranchLengthsForAllPartitions(MTInstance &instance) {
#if 0
    const unsigned numPartBranchLengths = instance.GetNumPartitions();
    instance.settings.partitionConverged.assign(numPartBranchLengths, 0);
    unsigned maxLoops = instance.settings.maxIterBrLenSmoothing;
    while (--maxLoops >= 0) {
        instance.settings.partitionSmoothed.assign(numPartBranchLengths, 1);
        sweepOverTreeOptimizeAllBranchLengthsForAllPartitions(instance);
        if (allSmoothed(tr, numPartBranchLengths))
            break;
    }
    instance.settings.partitionConverged.assign(numPartBranchLengths, 0);
#endif
}
#if 0
template<typename T>
class BeforeAfterIter {
    BeforeAfterIter(Node *, (*before)(Node *, T), (*after)(Node *, T), T blob);
};

class BeforeAfterIter {
    BeforeAfterIter(Node *, (*before)(Node *, void *), (*after)(Node *, void *), void * blob);
};

/**  sweepOverTreeOptimizeAllBranchLengthsForAllPartitions = smooth in PLL
*/
void sweepOverTreeOptimizeAllBranchLengthsForAllPartitions (MTInstance &instance)
{
    const unsigned numPartBranchLengths = instance.GetNumPartitions();
    Tree & tree = instance.tree;
    Node * root = tree.GetRoot();
    void * blob = &(instance);
    BeforeAfterIter edgeIter(root, optimizeAllLengthsForOneEdgeHook, updatCLAForExit, blob);
    while (edgeIter.get() != nullptr) {
        edgeIter.advance();
    }
}

boolean updatCLAForExit(Node * node, void * mtInstance) {
    if(numBranches > 1 && !tr->useRecom)
      pllUpdatePartials(tr, pr,p, PLL_TRUE);
    else
      pllUpdatePartials(tr, pr,p, PLL_FALSE);
}

//pll optimizeAllLengthsForOneEdge update
boolean optimizeAllLengthsForOneEdgeHook(Node * node, void * mtInstance) {
    MTInstance * instance = (MTInstance *)mtInstance;
    return optimizeAllLengthsForOneEdge(node, *instance);
}

boolean optimizeAllLengthsForOneEdgeHook(Arc edge, MTInstance & mtInstance) {
    const unsigned numPartBranchLengths = instance.GetNumPartitions();
    const double brLenDiffThreshold = mtInstance.settings.brLenDiffThreshold;
    vector<double> origEdgeLen(edge.GetEdgeLengthValue());
    vector<double> edgeLenCopy(origEdgeLen);
    makenewzGeneric(tr, p, q, z0, newzpercycle, &edgeLenCopy[0]);
    vector<char> smoothed = instance.settings.partitionSmoothed;
    for(auto i = 0U; i < numPartBranchLengths; i++) {
        if (instance.settings.partitionConverged[i]) {
            if(abs(edgeLenCopy[i] - origEdgeLen[i]) > brLenDiffThreshold) {
                smoothed[i] = 0;
            }
            edgeLengthPtr[i] = edgeLenCopy[i];
        }
    }
    instance.settings.partitionSmoothed = smoothed;
  return TRUE;
}


void update(pllInstance *tr, partitionList *pr, nodeptr p)
{
  nodeptr  q;
  int i;
  double   z[PLL_NUM_BRANCHES], z0[PLL_NUM_BRANCHES];
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;


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

/* the function below actually implements the iterative Newton-Raphson procedure.
   It is particularly messy and hard to read because for the case of per-partition branch length
   estimates it needs to keep track of whether the Newton Raphson procedure has
   converged for each partition individually.

   The rational efor doing it like this is also provided in:


   A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009,

*/

static void topLevelMakenewz(pllInstance *tr, partitionList * pr, double *z0, int _maxiter, double *result)
{
  double   z[PLL_NUM_BRANCHES], zprev[PLL_NUM_BRANCHES], zstep[PLL_NUM_BRANCHES];
  volatile double  dlnLdlz[PLL_NUM_BRANCHES], d2lnLdlz2[PLL_NUM_BRANCHES];
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
    tr->curvatOK[i]       = PLL_TRUE;
  }


  /* nested do while loops of Newton-Raphson */
  do
  {
    /* check if we ar done for partition i or if we need to adapt the branch length again */
    for(i = 0; i < numBranches; i++)
    {
      if(outerConverged[i] == PLL_FALSE && tr->curvatOK[i] == PLL_TRUE)
      {
        tr->curvatOK[i] = PLL_FALSE;

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
        lz    = log(z[i]);

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

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

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
#else
    /* sequential part, if this is the first newton-raphson implementation,
       do the precomputations as well, otherwise just execute the computation
       of the derivatives */
    if(firstIteration)
      {
        makenewzIterative(tr, pr);
        firstIteration = PLL_FALSE;
      }
    execCore(tr, pr, dlnLdlz, d2lnLdlz2);
    printf(" makenewz after execCore dlnLdlz = %lf d2lnLdlz2 = %lf\n", dlnLdlz, d2lnLdlz2);
#endif

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
    for(i = 0; i < numBranches; i++)
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

#endif

}// namespace
