#if ! defined(MT_OPT_MODEL_H)
#define MT_OPT_MODEL_H

namespace mt {

  class MTInstance;
  class CharModel;

  // definitely don't need all these functions in header file

  double LnGamma(double alph);

  double IncompleteGamma(double x, double alph, double lga);

  double PointNormal(double prob);

  double PointChi2(double prob, double df);

  void optimizeModel(MTInstance &instance, double likelihoodEpsilon);

  typedef std::pair<double,double> val_lnl_t;

  void optimizeModelUsingGolden(MTInstance &instance);

  static void changeParam(MTInstance &instance,
                          int ,//model,
                          int , //position,
                          double value,
                          int paramType);

  static void mtEvaluateChange(MTInstance &instance, int rateNumber, double value, double &result,
                               bool ,//converged,
                               int paramType, int model);

  static void mtBrentGeneric (MTInstance &instance, double &ax, double &bx, double &cx, double &fb, double tol,
                              double &xmin, double &result, int paramType, int rateNumber, double lowerBound,
                              double upperBound, int model);

  static int mtBrakGeneric(MTInstance &instance, double &param, double &ax, double &bx, double &cx, double &fa, double &fb,
                           double &fc, double lowerBound, double upperBound, int model, int rateNumber, int paramType);

  static void optParam(MTInstance &instance, int numberOfModels,
                       int rateNumber, double lowerBound, double upperBound, int paramType, double modelEpsilon);

  static void mtOptAlphas(MTInstance &instance, int numModels, double modelEpsilon);

  void initParams(MTInstance &instance);

  double fullOptimize(MTInstance &instance);

} // namespace mt

#endif
