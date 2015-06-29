#include "mt_char_model.h"
#include "mt_instance.h"
#include "mt_data.h"
#include "mt_tree.h"
#include "mt_tree_traversal.h"
#include <vector>

namespace mt {

// Parameter Optimization functions based on algorithms from "optimizeModel.c" in PLL, currently implemented for rate optimization, can be generalized for other params

// Parameter Type Values
#define SUB_RATE 0
#define ALPHA_P 1

// based on "changeModelParameters" - change param of type paramType to value at position
static void changeParam (MTInstance &instance, int position, double value, int paramType) {
  switch (paramType){

    case SUB_RATE:
      instance.changeRate(position, value);
      break;

    //case ALPHA_P:
      // write function to change alpha (shape) parameter, if necessary
    //  break;
  }

}

// Evaluate change to a model in terms of likelihood
// evaluateChange in PLL
static void mtEvaluateChange(MTInstance &instance, int rateNumber, std::vector<double> &values, std::vector<double> &results,
                           bool converged, int paramType, int nModels)
{
  int pos;

  switch (paramType)
  {
    case SUB_RATE:
      for(pos = 0; pos < nModels; pos++)
      {
        assert(rateNumber != -1);
        changeParam(instance, rateNumber, values[pos], paramType);
        CharModel &chrmodel = instance.GetCharModel();
        results[pos] = ScoreTree(instance.partMat, instance.tree, chrmodel);
        break;
      }
    //case ALPHA_P:
    //  break;
  }
}

// Implementation of Brent's Method for One-Dimensional Parameter Optimization, based on brentGeneric
// Currently assumes one partition, can be altered to account for more

static void mtBrentGeneric (MTInstance &instance, std::vector<double> &ax, std::vector<double> &bx, std::vector<double> &cx,
                            std::vector<double> &fb, double tol, std::vector<double> &xmin, std::vector<double> &result,
                            int paramType, int rateNumber, double lowerBound, double upperBound, int numberOfModels)
{

  //initialize variables
  bool allConverged = false;
  int i, iter;
  std::vector<bool> converged;
  std::vector<double> e, d, a, b, x, w, v, fw, fv, fx, xm, tol1, tol2, fu, r, q, p, u, etemp;

  //initialize Values + check
  for(i = 0; i < numberOfModels; i++)
    {
      converged[i] = false;

      e[i] = d[i] = 0.0;

      a[i] = std::min(cx[i],ax[i]),
      b[i] = std::max(cx[i],ax[i]),
      x[i] = bx[i], w[i] = bx[i], v[i] = bx[i],
      fw[i] = fb[i], fv[i] = fb[i], fx[i] = fb[i];

      assert(a[i] >= lowerBound && a[i] <= upperBound);
      assert(b[i] >= lowerBound && b[i] <= upperBound);
      assert(x[i] >= lowerBound && x[i] <= upperBound);
      assert(v[i] >= lowerBound && v[i] <= upperBound);
      assert(w[i] >= lowerBound && w[i] <= upperBound);
    }

  // core iteration
  for (iter = 0; iter < MAX_ITERS; iter++) {
    allConverged = true;

    for(i = 0; i < numberOfModels && allConverged; i++)
      allConverged = allConverged && converged[i];
    if (allConverged)
      return;

    for(i = 0; i < numberOfModels; i++)
      {
        if(!converged[i])
        {
          assert(a[i] >= lowerBound && a[i] <= upperBound);
          assert(b[i] >= lowerBound && b[i] <= upperBound);
          assert(x[i] >= lowerBound && x[i] <= upperBound);
          assert(v[i] >= lowerBound && v[i] <= upperBound);
          assert(w[i] >= lowerBound && w[i] <= upperBound);

          xm[i] = 0.5 * (a[i] + b[i]);
          tol2[i] = 2.0 * (tol1[i] = tol * fabs(x[i]) + BZ_EPSILON);

          if (fabs(x[i] - xm[i]) <= (tol2[i] - 0.5 * (b[i] - a[i])))
            {
              result[i] = -fx[i];
              xmin[i] = x[i];
              converged[i] = true;
            }
          else
            {
              if (fabs(e[i]) > tol1[i])
                {
                  r[i] = (x[i] - w[i]) * (fx[i] - fv[i]);
                  q[i] = (w[i] - v[i]) * (fx[i] - fw[i]);
                  p[i] = (x[i] - v[i]) * q[i] - (x[i] - w[i]) * r[i];
                  q[i] = 2.0 * (q[i] - r[i]);
                  if (q[i] > 0.0)
                    p[i] = -p[i];
                  q[i] = fabs(q[i]);
                  etemp[i] = e[i];
                  e[i] = d[i];
                  if((fabs(p[i]) >= fabs(0.5 * q[i] * etemp[i]))
                  || (p[i] <= q[i] * (a[i]-x[i])) || (p[i] >= q[i] * (b[i] - x[i])))
                    d[i] = BRENT_VAR * (e[i] = (x[i] >= xm[i] ? a[i] - x[i] : b[i] - x[i]));
                  else
                    {
                      d[i] = p[i] / q[i];
                      u[i] = x[i] + d[i];
                      if (u[i] - a[i] < tol2[i] || b[i] - u[i] < tol2[i])
                        d[i] = (xm[i] - x[i] > 0.0 ? fabs(tol1[i]) : -fabs(tol1[i]));

                    }
                }
              else
                {
                  d[i] = BRENT_VAR * (e[i] = (x[i] >= xm[i] ? a[i] - x[i] : b[i] - x[i]));
                }
              u[i] = ((fabs(d[i]) >= tol1[i]) ? (x[i] + d[i]) : (x[i] + (tol1[i] > 0.0 ? fabs(d[i]) : -fabs(d[i]))));
            }
          if (!converged[i])
            assert(u[i] >= lowerBound && u[i] <= upperBound);
        }

    }

  mtEvaluateChange(instance, rateNumber, u, fu, allConverged, paramType, numberOfModels);

  for(i = 0; i < numberOfModels; i++)
  {
    if(!converged[i])
      {
        if(fu[i] <= fx[i])
          {
            if(u[i] >= x[i])
              a[i] = x[i];
            else
              b[i] = x[i];

          v[i] = w[i]; w[i] = x[i]; x[i] = u[i];
          fv[i] = fw[i]; fw[i] = fx[i]; fx[i] = fu[i];

          }
        else
          {
            if (u[i] < x[i])
              a[i] = u[i];
            else
              b[i] = u[i];

            if (fu[i] <= fw[i] || w[i] == x[i])
              {
                v[i] = w[i];
                w[i] = u[i];
                fv[i] = fw[i];
                fw[i] = fu[i];
              }
            else
              {
                if(fu[i] <= fv[i] || v[i] == x[i] || v[i] == w[i])
                  {
                    v[i] = u[i];
                    fv[i] = fu[i];
                  }
              }
          }

        assert(a[i] >= lowerBound && a[i] <= upperBound);
        assert(b[i] >= lowerBound && b[i] <= upperBound);
        assert(x[i] >= lowerBound && x[i] <= upperBound);
        assert(v[i] >= lowerBound && v[i] <= upperBound);
        assert(w[i] >= lowerBound && w[i] <= upperBound);
        assert(u[i] >= lowerBound && u[i] <= upperBound);
      }
    }
  }
}

// Bracketing Function for Brent's Algorithm
// brakGeneric in PLL
static int mtBrakGeneric(MTInstance &instance, std::vector<double> &param, std::vector<double> &ax, std::vector<double> &bx,
                        std::vector<double>& cx, std::vector<double>& fa, std::vector<double>& fb, std::vector<double>& fc,
                        double lowerBound, double upperBound, int numberOfModels, int rateNumber, int paramType)
{
  std::vector<double> ulim, r, q, fu, dum, temp, u;
  int i;
  std::vector<int> states, endStates;
  std::vector<bool> converged;
  bool allConverged = false;

  for(i = 0; i < numberOfModels; i++)
    converged[i] = false;

  for(i = 0; i < numberOfModels; i++)
    {
      states[i] = 0;
      endStates[i] = 0;

      u[i] = 0.0;

      param[i] = ax[i];

      if(param[i] > upperBound)
        param[i] = ax[i] = upperBound;

      if(param[i] < lowerBound)
        param[i] = ax[i] = lowerBound;

      assert(param[i] >= lowerBound && param[i] <= upperBound);
    }

      mtEvaluateChange(instance, rateNumber, param, fa, allConverged, paramType, numberOfModels);

  for(i = 0; i < numberOfModels; i++)
    {
      param[i] = bx[i];
      if(param[i] > upperBound)
        param[i] = bx[i] = upperBound;
      if(param[i] < lowerBound)
        param[i] = bx[i] = lowerBound;

      assert(param[i] >= lowerBound && param[i] <= upperBound);
    }

      mtEvaluateChange(instance, rateNumber, param, fb, allConverged, paramType, numberOfModels);

  for(i = 0; i < numberOfModels; i++)
    {
      if (fb[i] > fa[i])
        {
          dum[i] = ax[i]; ax[i] = bx[i]; bx[i] = dum[i];
          dum[i] = fa[i]; fa[i] = fb[i]; fb[i] = dum[i];
        }

      cx[i] = bx[i] + GOLDEN_RAT * (bx[i] - ax[i]);

      param[i] = cx[i];

      if(param[i] > upperBound)
        param[i] = cx[i] = upperBound;
      if(param[i] < lowerBound)
        param[i] = cx[i] = lowerBound;

      assert(param[i] >= lowerBound && param[i] <= upperBound);
    }


      mtEvaluateChange(instance, rateNumber, param, fc, allConverged, paramType, numberOfModels);

   while(1)
     {
       allConverged = true;

       for(i = 0; i < numberOfModels && allConverged; i++)
         allConverged = allConverged && converged[i];

         if(allConverged)
          {
            for(i = 0; i < numberOfModels; i++)
              {
                if(ax[i] > upperBound)
                  ax[i] = upperBound;
                if(ax[i] < lowerBound)
                  ax[i] = lowerBound;

                if(bx[i] > upperBound)
                  bx[i] = upperBound;
                if(bx[i] < lowerBound)
                  bx[i] = lowerBound;

                if(cx[i] > upperBound)
                  cx[i] = upperBound;
                if(cx[i] < lowerBound)
                  cx[i] = lowerBound;
              }

            return 0;

          }

       for(i = 0; i < numberOfModels; i++)
         {
           if(!converged[i])
             {
               switch(states[i])
                 {
                 case 0:
                   endStates[i] = 0;
                   if(!(fb[i] > fc[i]))
                     converged[i] = true;
                   else
                     {

                       if(ax[i] > upperBound)
                         ax[i] = upperBound;
                       if(ax[i] < lowerBound)
                         ax[i] = lowerBound;
                       if(bx[i] > upperBound)
                         bx[i] = upperBound;
                       if(bx[i] < lowerBound)
                         bx[i] = lowerBound;
                       if(cx[i] > upperBound)
                         cx[i] = upperBound;
                       if(cx[i] < lowerBound)
                         cx[i] = lowerBound;

                       r[i]=(bx[i]-ax[i])*(fb[i]-fc[i]);
                       q[i]=(bx[i]-cx[i])*(fb[i]-fa[i]);
                       u[i]=(bx[i])-((bx[i]-cx[i])*q[i]-(bx[i]-ax[i])*r[i])/
                         (2.0 * MTREE_SIGN(std::max(fabs(q[i]-r[i]),BRAK_VAR1),q[i]-r[i]));

                       ulim[i]=(bx[i])+BRAK_VAR2*(cx[i]-bx[i]);

                       if(u[i] > upperBound)
                         u[i] = upperBound;
                       if(u[i] < lowerBound)
                         u[i] = lowerBound;
                       if(ulim[i] > upperBound)
                         ulim[i] = upperBound;
                       if(ulim[i] < lowerBound)
                         ulim[i] = lowerBound;

                       if ((bx[i]-u[i])*(u[i]-cx[i]) > 0.0)
                         {
                           param[i] = u[i];
                           if(param[i] > upperBound)
                             param[i] = u[i] = upperBound;
                           if(param[i] < lowerBound)
                             param[i] = u[i] = lowerBound;
                           endStates[i] = 1;
                         }
                       else
                         {
                           if ((cx[i]-u[i])*(u[i]-ulim[i]) > 0.0)
                             {
                               param[i] = u[i];
                               if(param[i] > upperBound)
                                 param[i] = u[i] = upperBound;
                               if(param[i] < lowerBound)
                                 param[i] = u[i] = lowerBound;
                               endStates[i] = 2;
                             }
                           else
                             {
                               if ((u[i]-ulim[i])*(ulim[i]-cx[i]) >= 0.0)
                                 {
                                   u[i] = ulim[i];
                                   param[i] = u[i];
                                   if(param[i] > upperBound)
                                     param[i] = u[i] = ulim[i] = upperBound;
                                   if(param[i] < lowerBound)
                                     param[i] = u[i] = ulim[i] = lowerBound;
                                   endStates[i] = 0;
                                 }
                               else
                                 {
                                   u[i]=(cx[i])+GOLDEN_RAT*(cx[i]-bx[i]);
                                   param[i] = u[i];
                                   endStates[i] = 0;
                                   if(param[i] > upperBound)
                                     param[i] = u[i] = upperBound;
                                   if(param[i] < lowerBound)
                                     param[i] = u[i] = lowerBound;
                                 }
                             }
                         }
                     }
                   break;
                 case 1:
                   endStates[i] = 0;
                   break;
                 case 2:
                   endStates[i] = 3;
                   break;
                 default:
                   assert(0);
                 }
               assert(param[i] >= lowerBound && param[i] <= upperBound);
             }
         }

          mtEvaluateChange(instance, rateNumber, param, temp, allConverged, paramType, numberOfModels);

       for(i = 0; i < numberOfModels; i++)
         {
           if(!converged[i])
             {
               switch(endStates[i])
                 {
                 case 0:
                   fu[i] = temp[i];
                   ax[i] = bx[i]; bx[i] = cx[i]; cx[i] = u[i];
                   fa[i] = fb[i]; fb[i] = fc[i]; fc[i] = fu[i];
                   states[i] = 0;
                   break;
                 case 1:
                   fu[i] = temp[i];
                   if (fu[i] < fc[i])
                     {
                       ax[i] = (bx[i]);
                       bx[i] = u[i];
                       fa[i] = (fb[i]);
                       fb[i] = fu[i];
                       converged[i] = true;
                     }
                   else
                     {
                       if (fu[i] > fb[i])
                         {
                           assert(u[i] >= lowerBound && u[i] <= upperBound);
                           cx[i] = u[i];
                           fc[i] = fu[i];
                           converged[i] = true;
                         }
                       else
                         {
                           u[i] = (cx[i])+GOLDEN_RAT*(cx[i]-bx[i]);
                           param[i] = u[i];
                           if(param[i] > upperBound)
                              param[i] = u[i] = upperBound;
                           if(param[i] < lowerBound)
                              param[i] = u[i] = lowerBound;
                           states[i] = 1;
                         }
                     }
                   break;
                 case 2:
                   fu[i] = temp[i];
                   if (fu[i] < fc[i])
                     {
                       bx[i] = cx[i]; cx[i] = u[i]; u[i] = cx[i] + GOLDEN_RAT*(cx[i]-bx[i]);
                       states[i] = 2;
                     }
                   else
                     {
                       states[i] = 0;
                       ax[i] = bx[i]; bx[i] = cx[i]; cx[i] = u[i];
                       fa[i] = fb[i]; fb[i] = fc[i]; fc[i] = fu[i];
                     }
                   break;
                 case 3:
                   fb[i] = fc[i]; fc[i] = fu[i]; fu[i] = temp[i];
                   ax[i] = bx[i]; bx[i] = cx[i]; cx[i] = u[i];
                   fa[i] = fb[i]; fb[i] = fc[i]; fc[i] = fu[i];
                   states[i] = 0;
                   break;
                 default:
                   assert(0);
                 }
             }
         }
    }

   return(0);
}

// Generic function for optimizing a given parameter
// from optParamGeneric in PLL
static void optParam(MTInstance &instance, int numberOfModels,
                     int rateNumber, double lowerBound, double upperBound, int paramType,
                     double modelEpsilon)
{
  int pos;

  std::vector<double> startRates, startWeights, startExps, startValues, startLHs,
         endLHs, _a, _b, _c, _fa, _fb, _fc, _param, _x;

  CharModel &startModel = instance.GetCharModel();

  for(pos = 0; pos < numberOfModels; pos++) {
    startLHs[pos] = ScoreTree(instance.partMat, instance.tree, startModel);

    switch (paramType)
    {
      case SUB_RATE:
        startValues[pos] = startModel.GetRate(rateNumber);
        break;
      case ALPHA_P:
      // function to return starting alpha parameter
        break;
    }

    _a[pos] = startValues[pos] + 0.1;
    _b[pos] = startValues[pos] - 0.1;

    if (_a[pos] < lowerBound)
      _a[pos] = lowerBound;
      if (_a[pos] > upperBound)
      _a[pos] = upperBound;
      if (_b[pos] < lowerBound)
      _b[pos] = lowerBound;
      if (_b[pos] > upperBound)
      _b[pos] = upperBound;


      mtBrakGeneric(instance, _param, _a, _b, _c, _fa, _fb, _fc, lowerBound, upperBound,
                  numberOfModels, rateNumber, paramType);

      assert(_a[pos] >= lowerBound && _a[pos] <= upperBound);
      assert(_b[pos] >= lowerBound && _b[pos] <= upperBound);
      assert(_c[pos] >= lowerBound && _c[pos] <= upperBound);

      mtBrentGeneric(instance, _a, _b, _c, _fb, modelEpsilon, _x, endLHs, paramType,
                 rateNumber, lowerBound, upperBound, numberOfModels);

      if (startLHs[pos] > endLHs[pos])
        {
          changeParam(instance, rateNumber, startValues[pos], paramType);
        }
      else {
          changeParam(instance, rateNumber, _x[pos], paramType);
        }
  }
}


//Wrapper function for optimizing all rates
// PLL's optRates
static void mtOptRates(MTInstance &instance, int numberOfModels, int states, double modelEpsilon)
{
  int
    rateNumber,
    numberOfRates = ((states * states - states) / 2) - 1;

  for(rateNumber = 0; rateNumber < numberOfRates; rateNumber++)
    optParam(instance, numberOfModels, rateNumber, MT_RATE_MIN, MT_RATE_MAX, SUB_RATE, modelEpsilon);
}


} //namespace mt
