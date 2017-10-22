/* This is a specialized implementation of the bootstrap to compute confidence intervals
   for angle deviations computed by dot products of suite of multiwavelet particle motion estimates.   

   Use this one for log amplitudes.
   */
class Bootstrap3DStatistics
{
  public:
    Bootstrap3DStatistics(double *x,int nx,double confidence, int number_trials);
    vector<double> median_vector()
    {
      return med;
    };
    double angle_confidence_interval()
    {
      return aci;
    };
    double confidence_level()
    {
      return cl;
    };
  private:
    /*! bootstrap median of 3D vectors */
    vector<double>  med;
    /* Estimated confidence interval.*/
    double aci;
    /* input confidence level */
    double cl;
};
/*! Generate random integers between 0 and nrange-1 */
int random_array_index(int range);
