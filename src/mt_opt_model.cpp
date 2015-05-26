#include "mt_char_model.h"
#include "mt_instance.h"
#include "mt_data.h"
#include "mt_tree.h"
#include "mt_tree_traversal.h"

using namespace std;
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
static void evaluateChange(MTInstance &instance, int rateNumber, double value, double result,
                           bool converged, int paramType, int numberOfModels, double modelEpsilon)
{
  switch (paramType)
  {
    case SUB_RATE:
      assert(rateNumber != -1);
      changeParam(instance, rateNumber, value, paramType);
      CharModel &chrmodel = instance.GetCharModel();
      result = ScoreTree(instance.partMat, instance.tree, chrmodel);
      break;
    //case ALPHA_P:
    //  break;
  }
}

// Implementation of Brent's Method for One-Dimensional Parameter Optimization, based on brentGeneric
// Currently assumes one partition, can be altered to account for more

static void mtBrentGeneric (MTInstance &instance, double ax, double bx, double cx, double fb, double tol,
                               double xmin, double result, int paramType, int position, double lowerBound,
                               double upperBound)
{

  //initialize variables
  bool converged = false;
  double e = 0.0,
         d = 0.0,
         a = std::min(cx,ax),
         b = std::max(cx,ax),
         x = bx, w = bx, v = bx,
         fw = fb, fv = fb, fx = fb;

  double xm, tol1, tol2, fu, r, q, p, u, etemp;

  assert(a >= lowerBound && a <= upperBound);
  assert(b >= lowerBound && b <= upperBound);
  assert(x >= lowerBound && x <= upperBound);
  assert(v >= lowerBound && v <= upperBound);
  assert(w >= lowerBound && w <= upperBound);

  // core iteration
  for (int i = 0; i < MAX_ITERS; i++) {

    if (converged)
      return;

    else {
      assert(a >= lowerBound && a <= upperBound);
      assert(b >= lowerBound && b <= upperBound);
      assert(x >= lowerBound && x <= upperBound);
      assert(v >= lowerBound && v <= upperBound);
      assert(w >= lowerBound && w <= upperBound);

      xm = 0.5 * (a + b);
      tol2 = 2.0 * (tol1 = tol * fabs(x) + BZ_EPSILON);

      if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
        {
          result = -fx;
          xmin = x;
          converged = true;
        }
      else
        {
          if (fabs(e) > tol1)
            {
              r = (x - w) * (fx - fv);
              q = (w - v) * (fx - fw);
              p = (x - v) * q - (x - w) * r;
              q = 2.0 * (q - r);
              if (q > 0.0)
                p = -p;
              q = fabs(q);
              etemp = e;
              e = d;
              if((fabs(p) >= fabs(0.5 * q * etemp))
                  || (p <= q * (a-x)) || (p >= q * (b - x)))
                d = BRENT_VAR * (e = (x >= xm ? a - x : b - x));
              else
                {
                  d = p / q;
                  u = x + d;
                  if (u - a < tol2 || b - u < tol2)
                    d = (xm - x > 0.0 ? fabs(tol1) : -fabs(tol1));

                }
            }
          else
            {
              d = BRENT_VAR * (e = (x >= xm ? a - x : b - x));
            }
          u = ((fabs(d) >= tol1) ? (x + d) : (x + (tol1 > 0.0 ? fabs(d) : -fabs(d))));
        }
    }
  }

  // write function for "evaluatechange"

  if(!converged)
    {
      if(fu <= fx)
        {
          if(u >= x)
            a = x;
          else
            b = x;

          v = w; w = x; x = u;
          fv = fw; fw = fx; fx = fu;

        }
      else
        {
          if (u < x)
            a = u;
          else
            b = u;

          if (fu <= fw || w == x)
            {
              v = w;
              w = u;
              fv = fw;
              fw = fu;
            }
          else
            {
              if(fu <= fv || v == x || v == w)
                {
                  v = u;
                  fv = fu;
                }
            }
        }

        assert(a >= lowerBound && a <= upperBound);
        assert(b >= lowerBound && b <= upperBound);
        assert(x >= lowerBound && x <= upperBound);
        assert(v >= lowerBound && v <= upperBound);
        assert(w >= lowerBound && w <= upperBound);
        assert(u >= lowerBound && u <= upperBound);
    }

}

// Bracketing Function for Brent's Algorithm
// brakGeneric in PLL
static int mtBrakGeneric(MTInstance &instance, double param, double ax, double bx, double cx, double fa, double fb,
                       double fc, double lowerBound, double upperBound,
                       int numberOfModels, int rateNumber, int whichFunction,
                       double modelEpsilon)
{
  double ulim, u, r, q, fu, dum, temp;

  int i, state, endState;

  bool converged = false;

  //for(i = 0; i < numberOfModels; i++)
  //  {
      state = 0;
      endState = 0;

      u = 0.0;

      param = ax;

      if(param > upperBound)
        param = ax = upperBound;

      if(param < lowerBound)
        param = ax = lowerBound;

      assert(param >= lowerBound && param <= upperBound);
  //  }

// to do: write evaluateChange function
//  evaluateChange(tr, pr, rateNumber, param, fa, converged, whichFunction, numberOfModels, ll, modelEpsilon);


//  for(i = 0; i < numberOfModels; i++)
//    {
      param = bx;
      if(param > upperBound)
        param = bx = upperBound;
      if(param < lowerBound)
        param = bx = lowerBound;

      assert(param >= lowerBound && param <= upperBound);
//    }

//  evaluateChange(tr, pr, rateNumber, param, fb, converged, whichFunction, numberOfModels, ll, modelEpsilon);

//  for(i = 0; i < numberOfModels; i++)
//    {
      if (fb > fa)
        {
          dum = ax; ax = bx; bx = dum;
          dum = fa; fa = fb; fb = dum;
        }

      cx = bx + GOLDEN_RAT * (bx - ax);

      param = cx;

      if(param > upperBound)
        param = cx = upperBound;
      if(param < lowerBound)
        param = cx = lowerBound;

      assert(param >= lowerBound && param <= upperBound);
  //  }


  // evaluateChange(tr, pr, rateNumber, param, fc, converged, whichFunction, numberOfModels,  ll, modelEpsilon);

   while(1)
     {

      // for(i = 0; i < numberOfModels && allConverged; i++)
    //     allConverged = allConverged && converged[i];

       if(converged)
         {
          // for(i = 0; i < numberOfModels; i++)
          //   {
               if(ax > upperBound)
                 ax = upperBound;
               if(ax < lowerBound)
                 ax = lowerBound;

               if(bx > upperBound)
                 bx = upperBound;
               if(bx < lowerBound)
                 bx = lowerBound;

               if(cx > upperBound)
                 cx = upperBound;
               if(cx < lowerBound)
                 cx = lowerBound;
          //   }

           return 0;

         }

      // for(i = 0; i < numberOfModels; i++)
      //   {
           if(!converged)
             {
               switch(state)
                 {
                 case 0:
                   endState = 0;
                   if(!(fb > fc))
                     converged = true;
                   else
                     {

                       if(ax > upperBound)
                         ax = upperBound;
                       if(ax < lowerBound)
                         ax = lowerBound;
                       if(bx > upperBound)
                         bx = upperBound;
                       if(bx < lowerBound)
                         bx = lowerBound;
                       if(cx > upperBound)
                         cx = upperBound;
                       if(cx < lowerBound)
                         cx = lowerBound;

                       r=(bx-ax)*(fb-fc);
                       q=(bx-cx)*(fb-fa);
                       u=(bx)-((bx-cx)*q-(bx-ax)*r)/
                         (2.0 * MTREE_SIGN(std::max(fabs(q-r),BRAK_VAR1),q-r));

                       ulim=(bx)+BRAK_VAR2*(cx-bx);

                       if(u > upperBound)
                         u = upperBound;
                       if(u < lowerBound)
                         u = lowerBound;
                       if(ulim > upperBound)
                         ulim = upperBound;
                       if(ulim < lowerBound)
                         ulim = lowerBound;

                       if ((bx-u)*(u-cx) > 0.0)
                         {
                           param = u;
                           if(param > upperBound)
                             param = u = upperBound;
                           if(param < lowerBound)
                             param = u = lowerBound;
                           endState = 1;
                         }
                       else
                         {
                           if ((cx-u)*(u-ulim) > 0.0)
                             {
                               param = u;
                               if(param > upperBound)
                                 param = u = upperBound;
                               if(param < lowerBound)
                                 param = u = lowerBound;
                               endState = 2;
                             }
                           else
                             {
                               if ((u-ulim)*(ulim-cx) >= 0.0)
                                 {
                                   u = ulim;
                                   param = u;
                                   if(param > upperBound)
                                     param = u = ulim = upperBound;
                                   if(param < lowerBound)
                                     param = u = ulim = lowerBound;
                                   endState = 0;
                                 }
                               else
                                 {
                                   u=(cx)+GOLDEN_RAT*(cx-bx);
                                   param = u;
                                   endState = 0;
                                   if(param > upperBound)
                                     param = u = upperBound;
                                   if(param < lowerBound)
                                     param = u = lowerBound;
                                 }
                             }
                         }
                     }
                   break;
                 case 1:
                   endState = 0;
                   break;
                 case 2:
                   endState = 3;
                   break;
                 default:
                   assert(0);
                 }
               assert(param >= lowerBound && param <= upperBound);
             }
      //   }

       // evaluateChange(tr, pr, rateNumber, param, temp, converged, whichFunction, numberOfModels, ll, modelEpsilon);

      // for(i = 0; i < numberOfModels; i++)
      //   {
           if(!converged)
             {
               switch(endState)
                 {
                 case 0:
                   fu = temp;
                   ax = bx; bx = cx; cx = u;
                   fa = fb; fb = fc; fc = fu;
                   state = 0;
                   break;
                 case 1:
                   fu = temp;
                   if (fu < fc)
                     {
                       ax=(bx);
                       bx=u;
                       fa=(fb);
                       fb=fu;
                       converged = true;
                     }
                   else
                     {
                       if (fu > fb)
                         {
                           assert(u >= lowerBound && u <= upperBound);
                           cx=u;
                           fc=fu;
                           converged = true;
                         }
                       else
                         {
                           u=(cx)+GOLDEN_RAT*(cx-bx);
                           param = u;
                           if(param > upperBound) {param = u = upperBound;}
                           if(param < lowerBound) {param = u = lowerBound;}
                           state = 1;
                         }
                     }
                   break;
                 case 2:
                   fu = temp;
                   if (fu < fc)
                     {
                       bx = cx; cx = u; u = cx + GOLDEN_RAT*(cx-bx);
                       state = 2;
                     }
                   else
                     {
                       state = 0;
                       ax = bx; bx = cx; cx = u;
                       fa = fb; fb = fc; fc = fu;
                     }
                   break;
                 case 3:
                   fb = fc; fc = fu; fu = temp;
                   ax = bx; bx = cx; cx = u;
                   fa = fb; fb = fc; fc = fu;
                   state = 0;
                   break;
                 default:
                   assert(0);
                 }
             }
      //   }
    }

   return(0);
}

// Generic function for optimizing a given parameter
// from optParamGeneric in PLL
static void optParam(MTInstance &instance, double modelEpsilon, int numberOfModels,
                     int rateNumber, double lowerBound, double upperBound, int paramType)
{
  int l,j,k,pos;

  double startRate, startWeight, startExp, startValue, startLH,
         endLH, _a, _b, _c, _fa, _fb, _fc, _param, _x;

  CharModel &startModel = instance.GetCharModel();
  startLH = ScoreTree(instance.partMat, instance.tree, startModel);

  switch (paramType)
  {
    case SUB_RATE:
      startValue = startModel.GetRate(rateNumber);
      break;
    case ALPHA_P:
      // function to return starting alpha parameter
      break;
  }

  _a = startValue + 0.1;
  _b = startValue - 0.1;

  if (_a < lowerBound)
    _a = lowerBound;
  if (_a > upperBound)
    _a = upperBound;
  if (_b < lowerBound)
    _b = lowerBound;
  if (_b > upperBound)
    _b = upperBound;


  mtBrakGeneric(instance, _param, _a, _b, _c, _fa, _fb, _fc, lowerBound, upperBound,
                numberOfModels, rateNumber, paramType, modelEpsilon);

  assert(_a >= lowerBound && _a <= upperBound);
  assert(_b >= lowerBound && _b <= upperBound);
  assert(_c >= lowerBound && _c <= upperBound);

  mtBrentGeneric(instance, _a, _b, _c, _fb, modelEpsilon, _x, endLH, paramType,
                 rateNumber, lowerBound, upperBound);

  if (startLH > endLH)
    {
      changeParam(instance, rateNumber, startValue, paramType);
    }
  else {
      changeParam(instance, rateNumber, _x, paramType);
    }
}

} //namespace mt
