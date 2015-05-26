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

    case ALPHA_P:
      // write function to change alpha (shape) parameter, if necessary
      break;
  }

}

// Implementation of Brent's Method for One-Dimensional Parameter Optimization, based on brentGeneric
// Currently assumes one partition, can be altered to account for more

static void mtreeBrentGeneric (MTInstance &instance, double ax, double bx, double cx, double fb, double tol,
                               double xmin, int paramType, int position, double lowerBound, double result,
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

} //namespace mt
