#include <vector>
#include "dmatrix.h"
/* This is a specialized implementation of the bootstrap to compute confidence intervals
   for angle deviations computed by dot products of suite of multiwavelet particle motion estimates.

   Use this one for log amplitudes.
   */
class Vector3DBootstrapError
{
  public:
    /*! Construct from a matrix of input values.

    This is the primary constructor for this object.  It uses construction
    is initialization, which in this case means the constructor does the
    bootstrap error estimation.  Estimates bootstrap errors at confidence
    level specified.  This program uses a delete and replace algorithm for
    every resampling.  i.e. randomly grab one vector from the trial set
    number_of_trials times with replacement after each draw.

    \param x - input data (assume sample vectors are in columns of x)
    \param - conficence is the confidence level computed. This must be a
    number greater than 0 and less than 1 or an error is thrown.
    \number_trials - number of resampling trials for the bootstrap.
    */
    Vector3DBootstrapError(dmatrix& x, const double confidence, const int number_trials);
    vector<double> mean_vector()
    {
      return mean;
    };
    double angle_error()
    {
      return aci;
    };
    double confidence_level()
    {
      return cl;
    };
  private:
    /*! bootstrap median of 3D vectors */
    vector<double>  mean;
    /* Estimated confidence interval.*/
    double aci;
    /* input confidence level */
    double cl;
};
/*! Generate random integers between 0 and nrange-1 */
int random_array_index(int range);
/* Simple procedure to estimate bootstrap mean and variance (the mv appendage)
for a vector of input numbers x.   ci is confidence level and ntrials is
number of trials.   Returns estiamte of center as first of pair and confidence
interval in second.  Note return is range divided by 2 to be comparable to a
sigma level */
pair<double,double> bootstrap_mv(const vector<double>& x, 
    const double ci, const double ntrials);
