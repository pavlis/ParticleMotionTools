#include <stdlib.h>
#include <stdio.h>
#include <vector>
/* these are needed to generate normally distributed random numbers */
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include "../Vector3DBootstrapError.h"
#include "dmatrix.h"
using namespace std;
using namespace boost;
double SampleNormal(const double mean, const double sigma)
{
    /* Generate a random number from the current time */
    //static mt19937 rng(static_cast<unsigned> (std::time(0)));
    uint64_t seed;
    seed=static_cast<uint64_t>(std::time(NULL));
    static boost::random::mt19937 rng(seed);
    /* This defines a normal distribution to generate numbers */
    boost::random::normal_distribution<double> norm_dist(mean,sigma);
    // bind random number generator to distribution to forma function object
    variate_generator<boost::random::mt19937&, boost::random::normal_distribution<double> >
        normal_sampler(rng, norm_dist);
    return(normal_sampler());
}
int main(int argc, char **argv)
{
  /* First create a set of normal deviates with mean 10 and sigma 1.  
     Use a most sample size */
  const int nx(100);
  const double xmean(10.0);
  const double xsigma(1.0);
  vector<double> x;
  x.reserve(nx);
  int i;
  cout << "Generating "<<nx<<" normal deviates with mean="<<xmean
    << " and sigma="<<xsigma<<endl;
  for(i=0;i<nx;++i)
  {
    x.push_back(SampleNormal(xmean,xsigma));
  }
  pair<double,double> bsresult;
  const double ci(0.95);
  const double ntrials(5000);
  bsresult=bootstrap_mv(x,ci,ntrials);
  cout << "Bootstrap mean="<<bsresult.first<<" 95% ci="<<bsresult.second<<endl;
  cout << "Now testing Vector3DBoostrapError class"<<endl;
  dmatrix d(3,nx);
  cout << "Generating "<<nx<<" 3D vectors with mean=(1,1,1) and "
    <<" sigma of each component=0.05"<<endl;
  for(i=0;i<nx;++i)
  {
    for(int k=0;k<3;++k)
      d(k,i)=SampleNormal(1.0,0.05);
  }
  cout << "Testing Vector3DBootstrapError constructor"<<endl;
  Vector3DBootstrapError dstat(d,ci,ntrials);
  cout << "Constructor finished"<<endl;
  vector<double> deltad;
  deltad=dstat.mean_vector();
  cout << "Bootstrap mean 3d vector=("<<deltad[0]<<", "<<deltad[1]
    <<", "<<deltad[2]<<")"<<endl;
  double dtheta=dstat.angle_error();
  cout << "Angle error (radians)="<<dtheta<<endl;
}
