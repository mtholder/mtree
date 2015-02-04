//Tree Search Algorithm from PLL - header file

#if !defined(__SEARCH_H__)
#define __SEARCH_H__
#include "mt_tree.h"

class rearrangeStep {
    public:
      enum rearrangeType {
          SPR_MODE = 0,
          NNI_MODE = 1
          // TBR_MODE = 2
      };
      double  likelihood;
      union   rearrangeType {
        class SPR {
          Node * removeNode;
          Node * insertNode;
          double zqr[NUM_BRANCHES];
          };
        class NNI {
            Node * originNode;
            int    swapType;
          };
      /*class TBR {
      }*/
      };
};

class nni_Move {
        MTInstance * tr;
        Node * p;
        int nniType;
        double z[PLL_NUM_BRANCHES]; // optimize branch lengths
        double z0[PLL_NUM_BRANCHES]; // unoptimized branch lengths
        double likelihood;
        double deltaLH;
};

class connectRELL {
  double z[PLL_NUM_BRANCHES];
  Node *p, *q;
  int cp, cq;
};

class topoRELL
{
  connectRELL     *connect;
  int             start;
  double          likelihood;
};


//Funcs

/* rearrange functions (NNI and SPR) */
topoList * createList (int max);
void topoSearch (MTInstance * tr, partitionList * pr, int rearrangeType, Node * p, int mintrav, int maxtrav, topoList * bestList);
void rearrangeCommit (MTInstance * tr, partitionList * pr, rearrangeInfo * rearr, int saveRollbackInfo);
int rearrangeRollback (MTInstance * tr, partitionList * pr);
void clearRearrangeHistory (MTInstance * tr);
int mtRaxmlSearchAlgorithm (MTInstance * tr, partitionList * pr, bool estimateModel);
int mtGetTransitionMatrix (MTInstance * tr, partitionList * pr, Node * p, int model, int rate, double * outBuffer);
void mtGetTransitionMatrix2 (MTInstance * tr, partitionList * pr, int model, Node * p, double * outBuffer);
int mtGetCLV (MTInstance * tr, partitionList * pr, Node * p, int partition, double * outProbs);
extern int mtTopologyPerformNNI(MTInstance * tr, Node * p, int swap);
