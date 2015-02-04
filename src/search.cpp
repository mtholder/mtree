#include "search.h"
#include "mt_tree.h"

#include <assert.h>

// Implementation of Tree Topology Search Algorithm ported from PLL

void performSearch (MTInstance &instance, int steps, Node *p);



static void saveTopolRELLRec(MTInstance *tr, Node *p, topolRELL *tpl, int *i, int numsp)
{
  int k;
  if(isTip(p->number, numsp))
    return;
  else
    {
      nodeptr q = p->next;
      while(q != p)
    {
      tpl->connect[*i].p = q;
      tpl->connect[*i].q = q->back;

      if(tr->grouped ||  tr->constrained)
        {
          tpl->connect[*i].cp = tr->constraintVector[q->number];
          tpl->connect[*i].cq = tr->constraintVector[q->back->number];
        }

      for(k = 0; k < PLL_NUM_BRANCHES; k++)
        tpl->connect[*i].z[k] = q->z[k];
      *i = *i + 1;

      saveTopolRELLRec(tr, q->back, tpl, i, numsp);
      q = q->next;
    }
    }
}

static void saveTopolRELL(pllInstance *tr, topolRELL *tpl)
{
  nodeptr p = tr->start;
  int k, i = 0;

  tpl->likelihood = tr->likelihood;
  tpl->start      = 1;

  tpl->connect[i].p = p;
  tpl->connect[i].q = p->back;

  if(tr->grouped ||  tr->constrained)
    {
      tpl->connect[i].cp = tr->constraintVector[p->number];
      tpl->connect[i].cq = tr->constraintVector[p->back->number];
    }

  for(k = 0; k < PLL_NUM_BRANCHES; k++)
    tpl->connect[i].z[k] = p->z[k];
  i++;

  saveTopolRELLRec(tr, p->back, tpl, &i, tr->mxtips);

  assert(i == 2 * tr->mxtips - 3);
}
