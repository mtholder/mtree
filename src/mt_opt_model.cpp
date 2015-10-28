#include "mt_char_model.h"
#include "mt_instance.h"
#include "mt_data.h"
#include "mt_opt_model.h"
#include "mt_tree.h"
#include "mt_tree_traversal.h"
#include "mt_likelihood.h"
#include <vector>
#include <cmath>

namespace mt {

// Parameter Optimization functions based on algorithms from "optimizeModel.c" in PLL, currently implemented for rate optimization, can be generalized for other params

// Parameter Type Values
#define SUB_RATE 0
#define ALPHA_P 1



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
  return f + (a-0.5)*log(a) - a + .918938533204673
         + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
         +.083333333333333)/a;
}

// Incomplete gamma ratio, x = upper limit of integration
// Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
// 19: 285-287 (AS32)
double IncompleteGamma(double x, double alph, double lga){
  int i;
  double p = alph, g = lga, accurate = 1e-8, overflow = 1e30, gin = 0, rn = 0;
  double a = 0, b = 0, an = 0, dif = 0, term = 0, pn[6], factor;

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
   for (i = 0; i < 4; i++) pn[i] /= overflow;
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

  assert(prob < 0.999998 && prob > 0.000002 && df > 0);

  double e = .5e-6,
         aa = .6931471805,
         p = prob,
         g, xx, c, ch,
         a = 0, q = 0, p1 = 0, p2 = 0, t = 0, x = 0, b = 0,
         s1, s2, s3, s4, s5, s6;

  g = LnGamma(df/2);
  //_DEBUG_VAL(g);
  xx = df/2;
  c = xx - 1;
  if(df >= -1.24*log(p)) goto l1;
  ch = pow((p*xx*exp(g+xx*aa)), 1/xx);
  //_DEBUG_VAL(ch);
  if (ch-e<0) return ch;
  goto l4;

l1:
  if (df > .32) goto l3;
  ch = 0.4;
  a = log(1 - p);
  //_DEBUG_VAL(a);

l2:
  q = ch;
  p1 = 1 + ch*(4.67 + ch);
  p2 = ch*(6.73 + ch*(6.66 + ch));
  t = -0.5 + (4.67 + 2*ch)/p1 - (6.73 + ch*(13.32 + 3*ch))/p2;
  ch -= (1 - exp(a + g + 0.5*ch + c*aa)*p2/p1)/t;
  //_DEBUG_VAL(ch);
  if(fabs(q/ch - 1) - 0.1 <= 0) goto l4;
  else goto l2;

l3:
  x = PointNormal(p);
  p1 = 0.222222/df;
  ch = df*pow((x*sqrt(p1) + 1 - p1), 3.0);
  if(ch > 2.2*df + 6)
    ch = -2*(log(1 - p) - c*log(0.5*ch) + g);

l4:
  q=ch;
  //_DEBUG_VAL(q);
  p1 = .5 * ch;
  t = IncompleteGamma(p1, xx, g);
  //_DEBUG_VAL(t);
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
  //_DEBUG_VAL(ch);
  if (fabs(q/ch - 1) > e) goto l4;

  //std::cerr << "Got here\n";
  return (ch);
}

// initialize general CharModel params
// Iterate along modelParams and initialize a set of parameters for each partition
/*void CharModel::initModels(PartitionedMatrix &partMat, unsigned modelType, std::vector<int> stateSetSizes) {
    for (int i = 0; i < partMat.GetNumPartitions(); i++) {
      ModelParams &modelParams(stateSetSizes[i], modelType);
      modelParams.initializeModel();
      this->modelList[i] = modelParams;
    }
}
*/

void CharModel::createGammas(double alph, int cats){
  int i;
  this->alpha = alph;
  //_DEBUG_VAL(alpha);
  double factor = alpha / alpha * cats;
  double beta = alpha;
  std::vector<double> gammaProbs(cats,0.0);

  assert(alpha >= MT_ALPHA_MIN);

  double lngal = LnGamma(alpha + 1);
  //_DEBUG_VAL(lngal);
  for(i = 0; i < cats - 1; i++) {
    gammaProbs[i] = MT_POINT_GAMMA((i + 1.0)/cats, alpha, beta);
    //_DEBUG_VAL(gammaProbs[i]);
  }
  for(i = 0; i < cats - 1; i++)
    gammaProbs[i] = IncompleteGamma(gammaProbs[i] * beta, alpha + 1, lngal);
  rates[0] = gammaProbs[0] * factor;
  rates[cats - 1] = (1 - gammaProbs[cats - 2]) * factor;
  for(i = 1; i < cats - 1; i++) {
    rates[i] = (gammaProbs[i] - gammaProbs[i - 1]) * factor;
    //_DEBUG_VAL(rates[i]);
  }
  return;
}

// based on "changeModelParameters" - change param of type paramType to value at position
// old value needs to be stored first!
static void changeParam(MTInstance &instance,
                        int model,
                        int , //position,
                        double value,
                        int paramType) {
  switch (paramType){
    case SUB_RATE:
      //instance.changeRate(model, position, value);
      //break;
    case ALPHA_P:
      instance.GetCharModel(model).createGammas(value, instance.GetCharModel(model).GetNumRates());
      //_DEBUG_VAL(value);
      break;
  }

}

// Evaluate change to a model in terms of likelihood
// equivalent to evaluateChange in PLL
static void mtEvaluateChange(MTInstance &instance, int rateNumber, double value, double &result,//result,
                           bool ,//converged,
                           int paramType, int model)
{
  switch (paramType)
  {
    case SUB_RATE:
      /*for(pos = 0; pos < nModels; pos++)
      {
        assert(rateNumber != -1);
        changeParam(instance, pos, rateNumber, values[pos], paramType);
        CharModel &chrmodel = instance.GetCharModel(0);
        results[pos] = ScoreTreeForPartition(instance.partMat, instance.tree, chrmodel, 0);
        break;
      }*/
    case ALPHA_P:
      changeParam(instance, model, rateNumber, value, paramType);
      break;
  }
  instance.likelihoods[model] = ScoreTreeForPartition(instance.partMat, instance.tree, instance.GetCharModel(model),model);
  result = instance.likelihoods[model];
  //_DEBUG_VAL(result);
  instance.dirtyFlags[model] = true;
}

// Implementation of Brent's Method for One-Dimensional Parameter Optimization, based on brentGeneric in PLL

void mtBrentGeneric (MTInstance &instance, double &ax, double &bx, double &cx, double &fb, double tol,
                            double &xmin, double &result, int paramType, int rateNumber, double lowerBound,
                            double upperBound, int model)
{

  //initialize variables
  bool allConverged = false;
  int iter;
  bool converged;
  double e, d, a, b, x, w, v, fw, fv, fx, xm, tol1, tol2, fu, r, q, p, u=0.0, etemp;

  //initialize Values + check
  converged = false;

  e = d = 0.0;

  a = std::min(cx,ax),
  b = std::max(cx,ax),
  x = bx, w = bx, v = bx,
  fw = fb, fv = fb, fx = fb;

  assert(a >= lowerBound && a <= upperBound);
  assert(b >= lowerBound && b <= upperBound);
  assert(x >= lowerBound && x <= upperBound);
  assert(v >= lowerBound && v <= upperBound);
  assert(w >= lowerBound && w <= upperBound);

  // core iteration
  for (iter = 0; iter < MAX_ITERS; iter++) {

    //std::cerr << "Brent iteration #: " << iter << "\n";
    if (converged)
      break; //MTH changed from return

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
      } else {
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
            || (p <= q * (a-x)) || (p >= q * (b - x))) {
              d = BRENT_VAR * (e = (x >= xm ? a - x : b - x));
            } else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                  d = (xm - x > 0.0 ? fabs(tol1) : -fabs(tol1));

              }
          } else {
              d = BRENT_VAR * (e = (x >= xm ? a - x : b - x));
            }
        u = ((fabs(d) >= tol1) ? (x + d) : (x + (tol1 > 0.0 ? fabs(d) : -fabs(d))));
      }
    if (!converged)
      assert(u >= lowerBound && u <= upperBound);
    //_DEBUG_VAL(xm);
    //std::cerr << "Current alpha: " << instance.GetCharModel(model).GetAlpha() << "\n";


  mtEvaluateChange(instance, rateNumber, u, fu, allConverged, paramType, model);


  if(!converged) {
     if(fu <= fx) {
       if(u >= x)
         a = x;
       else
         b = x;
         v = w; w = x; x = u;
         fv = fw; fw = fx; fx = fu;
      } else {
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
           } else {
             if(fu <= fv || v == x || v == w) {
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

        xmin = x;
      }
  }
}

// Bracketing Function for Brent's Algorithm
// brakGeneric in PLL
int mtBrakGeneric(MTInstance &instance, double &param, double &ax, double &bx, double &cx, double &fa, double &fb,
                        double &fc, double lowerBound, double upperBound, int model, int rateNumber, int paramType)
{
  //_DEBUG_VAL(cx);
  double ulim, r, q, fu=0.0, dum, temp, u;
  int state, endState;
  bool converged = false;

  state = 0;
  endState = 0;

  u = 0.0;

  param = ax;

  if(param > upperBound)
    param = ax = upperBound;

  if(param < lowerBound)
    param = ax = lowerBound;

  assert(param >= lowerBound && param <= upperBound);

  mtEvaluateChange(instance, rateNumber, param, fa, converged, paramType, model);

   param = bx;
   if(param > upperBound)
     param = bx = upperBound;
   if(param < lowerBound)
     param = bx = lowerBound;

   assert(param >= lowerBound && param <= upperBound);

   mtEvaluateChange(instance, rateNumber, param, fb, converged, paramType, model);

   if (fb > fa) {
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

   mtEvaluateChange(instance, rateNumber, param, fc, converged, paramType, model);

   while(1)
     {
         if(converged) {

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

            return 0;
          }

         if(!converged) {
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
                  } else {
                    if ((cx-u)*(u-ulim) > 0.0)
                    {
                      param = u;
                      if(param > upperBound)
                        param = u = upperBound;
                      if(param < lowerBound)
                        param = u = lowerBound;
                      endState = 2;
                    } else {
                      if ((u-ulim)*(ulim-cx) >= 0.0)
                      {
                        u = ulim;
                        param = u;
                        if(param > upperBound)
                          param = u = ulim = upperBound;
                        if(param < lowerBound)
                          param = u = ulim = lowerBound;
                        endState = 0;
                       } else {
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
                   _DEBUG_VAL(cx);
                   break;
                 case 1:
                   endState = 0;
                   _DEBUG_VAL(cx);
                   break;
                 case 2:
                   endState = 3;
                   _DEBUG_VAL(cx);
                   break;
                 default:
                   assert(0);
                 }
               assert(param >= lowerBound && param <= upperBound);
             }

          mtEvaluateChange(instance, rateNumber, param, temp, converged, paramType, model);

           if(!converged)
             {
               switch(endState)
                 {
                 case 0:
                   fu = temp;
                   ax = bx; bx = cx; cx = u;
                   fa = fb; fb = fc; fc = fu;
                   state = 0;
                   //_DEBUG_VAL(cx);
                   break;
                 case 1:
                   fu = temp;
                   if (fu < fc)
                     {
                       ax = (bx);
                       bx = u;
                       fa = (fb);
                       fb = fu;
                       converged = true;
                     }
                   else
                     {
                       if (fu > fb)
                         {
                           assert(u >= lowerBound && u <= upperBound);
                           cx = u;
                           fc = fu;
                           converged = true;
                         }
                       else
                         {
                           u = (cx)+GOLDEN_RAT*(cx-bx);
                           param = u;
                           if(param > upperBound)
                              param = u = upperBound;
                           if(param < lowerBound)
                              param = u = lowerBound;
                           state = 1;
                         }
                     }
                   //_DEBUG_VAL(cx);
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
                   //_DEBUG_VAL(cx);
                   break;
                 case 3:
                   fb = fc; fc = fu; fu = temp;
                   ax = bx; bx = cx; cx = u;
                   fa = fb; fb = fc; fc = fu;
                   state = 0;
                   //_DEBUG_VAL(cx);
                   break;
                 default:
                   assert(0);
                 }
             }
    }
//_DEBUG_VAL(cx);
}

// Set beginning parameters for optimization
void initParams(MTInstance &instance) {
  for(auto i = 0U; i < instance.GetModelVec().size(); i++) { // loop over partitions
    instance.GetCharModel(i).createGammas(0.5, instance.GetCharModel(i).GetNumRates());
    std::cerr << "Set alpha for partiton " << i << "\n";
    // Set alphas to 0.5 for now
  }
}


// Generic function for optimizing a given parameter
// from optParamGeneric in PLL
// loop over parameters occurs here - specific to parameter in nested functions
static void optParam(MTInstance &instance, int numberOfModels,
                     int rateNumber, double lowerBound, double upperBound, int paramType,
                     double modelEpsilon)
{
  int pos;

  std::vector<double> startRates, startWeights, startExps;

  for(pos = 0; pos < numberOfModels; pos++) {
    CharModel &startModel = instance.GetCharModel(pos);
    double startLH = ScoreTreeForPartition(instance.partMat, instance.tree, startModel, pos);
    //_DEBUG_VAL(startLH);
    double startValue = 0.0;
    switch (paramType)
    {
      case SUB_RATE:
        //startValues[pos] = startModel.GetRate(rateNumber);
        //break;
      case ALPHA_P:
        // function to return starting alpha parameter
        startValue = GetPatData(pos).GetAlpha();
        _DEBUG_VAL(startValue);
        break;
    }

    double _a = startValue + 0.1;
    double _b = startValue - 0.1;

    if (_a < lowerBound)
      _a = lowerBound;
      if (_a > upperBound)
      _a = upperBound;
      if (_b < lowerBound)
      _b = lowerBound;
      if (_b > upperBound)
      _b = upperBound;

      double _c, _fa, _fb, _fc, _param, _x, endLH;
      //_x = startValue;
      //std::cerr << "Got to mtBrakGeneric\n";
      mtBrakGeneric(instance, _param, _a, _b, _c, _fa, _fb, _fc, lowerBound, upperBound,
                    pos, rateNumber, paramType);

      //_DEBUG_VAL(_c);

      assert(_a >= lowerBound && _a <= upperBound);
      assert(_b >= lowerBound && _b <= upperBound);
      assert(_c >= lowerBound && _c <= upperBound);

      //_DEBUG_VAL(_x);
      mtBrentGeneric(instance, _a, _b, _c, _fb, modelEpsilon, _x, endLH, paramType,
                     rateNumber, lowerBound, upperBound, pos);

      //std::cerr << "Got through Brent\n";
      if (startLH > endLH)
        {
          changeParam(instance, pos, rateNumber, startValue, paramType);

        }
      else {
          //_DEBUG_VAL(_x);
          changeParam(instance, pos, rateNumber, _x, paramType);
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

static void mtOptAlphas(MTInstance &instance, int numModels, double modelEpsilon)
{
  //std::cerr << "Entering mtOptAlphas\n";
  optParam(instance, numModels, 0, 0.01, 0.99, ALPHA_P, modelEpsilon);
}

// iterative procedure for optimizing all model parameters for all partitions
// "modOpt" in optimizeModel.c in PLL
// for now only optimizes branch lengths, alpha params, and rates
void optimizeModel(MTInstance &instance, double likelihoodEpsilon) {
  std::cerr << "Beginning optimization\n";
  initParams(instance);
  double inputLikelihood, currentLikelihood, modelEpsilon;

  //std::vector<double> alphalist, ratelist; // get these from instance, where they have been initialized

  modelEpsilon = 0.0001;    // hard-coded value for now

  inputLikelihood = instance.curLikelihood;
  int MAX_ITS = 10;
  do {
    std::cerr << "Optimization iterations left: " << MAX_ITS << '\n';
    currentLikelihood = instance.curLikelihood;

    mtOptAlphas(instance, instance.numPartitions, modelEpsilon);
    instance.curLikelihood = ScoreTree(instance.partMat, instance.tree, instance);
    MAX_ITS--;

  } while (fabs(currentLikelihood - instance.curLikelihood) > likelihoodEpsilon && MAX_ITS);

}

} //namespace mt
