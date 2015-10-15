#include "mt_char_model.h"
#include "mt_instance.h"
#include "mt_data.h"
#include "mt_tree.h"
#include "mt_tree_traversal.h"
#include <vector>
#include <cmath>

namespace mt {

// Parameter Optimization functions based on algorithms from "optimizeModel.c" in PLL, currently implemented for rate optimization, can be generalized for other params

// Parameter Type Values
#define SUB_RATE 0
#define ALPHA_P 1

// Model Type Values
#define _MK_           0
#define _MK_COVAR_     1


// The following four functions(LnGamma, IncompleteGamma, PointNormal, and PointChi2)
// come from Ziheng Yang's PAML

// ln(gamma(alpha)) using Stirling's formula
// Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
// Communications of the Association for Computing Machinery, 9:684
double LnGamma(double alph){

  double a = alph,
         f = 0,
         z;

  if (a < 7) {
    f = 1;
    z = a - 1;
    while(++z < 7)
      f *= z;
    a = z;
    f = -log(f);
  }
  z = 1 / (a * a);
  return f + (a-0.5)*log(x) - a + .918938533204673
	       + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	       +.083333333333333)/a;
}

// Incomplete gamma ratio, x = upper limit of integration
// Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
// 19: 285-287 (AS32)
double IncompleteGamma(double x, double alph, double lga){
  int i;
  double p = alpha, g = lga, accurate = 1e-8, overflow = 1e30;
  double factor, gin = 0, rn = 0, a = 0, b = 0, an = 0, dif = 0, term = 0, pn[6];

  if(x == 0)
    return 0;
  assert(x >= 0 && p > 0);
  factor = exp(p*log(x) - x - g);
  if (x > 1 && x >= p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term *= x/rn;   gin += term;
   if (term > accurate) goto l20;
   gin *= factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a = 1 - p;   b = a + x + 1;  term = 0;
   pn[0] = 1;  pn[1] = x;  pn[2] = x + 1;  pn[3] = x*b;
   gin = pn[2]/pn[3];
 l32:
   a++;  b += 2;  term++;   an = a*term;
   for (i = 0; i < 2; i++) pn[i + 4]=b*pn[i + 2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn = pn[4]/pn[5];   dif = fabs(gin-rn);
   if (dif > accurate) goto l34;
   if (dif <= accurate*rn) goto l42;
 l34:
   gin = rn;
 l35:
   for (i = 0; i < 4; i++) pn[i] = pn[i + 2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i = 0; i < 4; i++) pn[i]/ = overflow;
   goto l32;
 l42:
   gin = 1 - factor*gin;

 l50:
   return (gin);
}

// Inverse distribution function for normal distribution
// Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
// Applied Statistics 22: 96-97 (AS70)
double PointNormal(double prob){
  double a0 = -.322232431088, a1 = -1, a2 = -.342242088547, a3 = -.0204231210245;
  double a4 = -.453642210148e-4, b0 = .0993484626060, b1 = .588581570495;
  double b2 = .531103462366, b3 = .103537752850, b4 = .0038560700634;
  double y, z = 0, p = prob, p1;

  p1 = (p<0.5 ? p : 1-p);
  assert(p1 > 1e-20);

  y = sqrt (log(1/(p1*p1)));
  z = y + ((((y*a4 + a3)*y + a2)*y + a1)*y + a0) / ((((y*b4 + b3)*y + b2)*y + b1)*y + b0);
  return (p < 0.5 ? -z : z);
}

// inverse distribution function of Chi distribution at df degrees of freedom
// returns z such that P{x < z} = prob
// note: uses goto + labels because imported from fortran code, alter for more readability
double PointChi2(double prob, double df){

  assert(prob < 0.999998 && prob > 0.000002 && v > 0);

  double e = .5e-6,
         aa = .6931471805,
         p = prob,
         g, xx, c, ch,
         a = 0, q = 0, p1 = 0, p2 = 0, t = 0, x = 0, b = 0,
         s1, s2, s3, s4, s5, s6;

  g = LnGamma(v/2);
  xx = v/2;
  c = xx - 1;
  if(v >= -1.24*log(p)) goto l1;
  ch = pow((p*xx*exp(g+xx*aa)), 1/xx);
  if (ch < e) return ch;
  goto l4;

l1:
  if (v > .32) goto l3;
  ch = 0.4;
  a = log(1 - p);

l2:
  q = ch;
  p1 = 1 + ch*(4.67 + ch);
  p2 = ch*(6.73 + ch*(6.66 + ch));
  t = -0.5 + (4.67 + 2*ch)/p1 - (6.73 + ch*(13.32 + 3*ch))/p2;
  ch -= (1 - exp(a + g + 0.5*ch + c*aa)*p2/p1)/t;
  if(fabs(q/ch - 1) - 0.1 <= 0) goto l4;
  else goto l2;

l3:
  x = PointNormal(p);
  p1 = 0.222222/v;
  ch = v*pow((x*sqrt(p1) + 1 - p1), 3.0);
  if(ch > 2.2*v + 6)
    ch = -2*(log(1 - p) - c*log(0.5*ch) + g);

l4:
  p1 = .5 * ch;
  t = IncompleteGamma(p1, xx, g);
  assert(t >= 0);
  p2 = p - t;
  t = p2*exp(xx*aa + g + p1 - c*log(ch));
  b = t/ch;
  a = .5*t - b*c;

  s1 = (210 + a*(140 + a*(105 + a*(84 + a*(70 + 60*a))))) / 420;
  s2 = (420 + a*(735 + a*(966 + a*(1141 + 1278*a))))/2520;
  s3 = (210 + a*(462 + a*(707 + 932*a)))/2520;
  s4 = (252 + a*(672 + 1182*a) + c*(294 + a*(889 + 1740*a)))/5040;
  s5 = (84 + 264*a + c*(175 + 606*a))/2520;
  s6 = (120 + c*(346 + 127*c))/5040;
  ch += t*(1 + 0.5*t*s1 - b*c*(s1 - b*(s2 - b*(s3 - b*(s4 - b*(s5 - b*s6))))));
  if (fabs(q/ch - 1) > e) goto l4;

  return (ch);
}

// initialize general CharModel params
// Iterate along modelParams and initialize a set of parameters for each partition
void CharModel::initModels(PartitionedMatrix &partMat, unsigned modelType, std::vector<int> stateSetSizes) {
    for (int i = 0; i < partMat.GetNumPartitions(); i++) {
      ModelParams &modelParams(stateSetSizes[i], modelType);
      modelParams.initializeModel();
      this->modelList[i] = modelParams;
    }
}

void ModelParams::createGammas(double alph, int cats){
  int i;
  alpha = alph;
  double factor = alpha / alpha * cats;
  double beta = alpha;
  std::vector<double> gammaProbs;

  assert(alpha >= MT_ALPHA_MIN);

  lngal = LnGamma(alpha + 1);
  for(i = 0; i < cat - 1; i++)
    gammaProbs[i] = MT_POINT_GAMMA((i + 1.0)/cat, alpha, beta);
  for(i = 0; i < cat - 1; i++)
    gammaProbs[i] = IncompleteGamma(gammaProbs[i] * beta, alpha + 1, lngal);
  gammaRates[0] = gammaProbs[0] * factor;
  gammaRates[cat - 1] = (1 - gammaProbs[cat - 2]) * factor;
  for(i = 1; i < cat - 1; i++)
    gammaRates[i] = (gammaProbs[i] - gammaProbs[i - 1]) * factor;

  return;
}

// based on "changeModelParameters" - change param of type paramType to value at position
// old value needs to be stored first!
static void changeParam(MTInstance &instance, int model, int position, double value, int paramType) {
  switch (paramType){

    case SUB_RATE:
      //instance.changeRate(model, position, value);
      //break;

    case ALPHA_P:
      instance.GetCharModel().GetModelParams(model).createGammas(value, instance.GetCharModel().GetNumRates());
      break;
  }

}

// Evaluate change to a model in terms of likelihood
// equivalent to evaluateChange in PLL
static void mtEvaluateChange(MTInstance &instance, int rateNumber, std::vector<double> values, std::vector<double> results,
                           bool converged, int paramType, int nModels)
{
  int pos;

  switch (paramType)
  {
    case SUB_RATE:
      for(pos = 0; pos < nModels; pos++)
      {
        assert(rateNumber != -1);
        changeParam(instance, pos, rateNumber, values[pos], paramType);
        CharModel &chrmodel = instance.GetCharModel();
        results[pos] = ScoreTree(instance.partMat, instance.tree, chrmodel);
        break;
      }
    case ALPHA_P:
      for(pos = 0; pos < nModels; pos++)
      {
      	changeParam(instance, pos, rateNumber, values[pos], paramType);
      	
      }
      break;
  }
}

// Implementation of Brent's Method for One-Dimensional Parameter Optimization, based on brentGeneric in PLL

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
    ModelParams &startParams = startModel.GetModelParams(pos);
    startLHs[pos] = ScoreTree(instance.partMat, instance.tree, startModel);

    switch (paramType)
    {
      case SUB_RATE:
        //startValues[pos] = startModel.GetRate(rateNumber);
        //break;
      case ALPHA_P:
      // function to return starting alpha parameter
      startValues[pos] = startParams.GetAlpha();
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

/*
//Wrapper function for optimizing all rates
// PLL's optRates
static void mtOptRates(MTInstance &instance, std::vector<double> ratelist, double modelEpsilon)
{
  int
    rateNumber,
    numberOfRates = ((states * states - states) / 2) - 1,
    numberOfModels = instance.partMat.GetNumPartitions();

  for(rateNumber = 0; rateNumber < numberOfRates; rateNumber++)
    optParam(instance, numberOfModels, rateNumber, MT_RATE_MIN, MT_RATE_MAX, SUB_RATE, modelEpsilon);
}
*/

static void mtOptAlphas(MTInstance &instance, std::vector<double> alphaList, double modelEpsilon)
{

}

// iterative procedure for optimizing all model parameters for all partitions
// "modOpt" in optimizeModel.c in PLL
// for now only optimizes branch lengths, alpha params, and rates
void optimizeModel(MTInstance &instance, double likelihoodEpsilon) {
  int i = 0;
  double inputLikelihood, currentLikelihood, modelEpsilon;

  std::vector<double> alphalist, ratelist; // get these from instance, where they have been initialized

  modelEpsilon = 0.0001;    // hard-coded value for now

  inputLikelihood = instance.curLikelihood;
  do {
    currentLikelihood = instance.curLikelihood;

    mtOptAlphas(instance, ratelist, modelEpsilon);
    instance.curLikelihood = ScoreTree(instance.partMat, instance.tree, instance.GetCharModel());

  } while (fabs(currentLikelihood - instance.curLikelihood) > likelihoodEpsilon)

}

} //namespace mt

