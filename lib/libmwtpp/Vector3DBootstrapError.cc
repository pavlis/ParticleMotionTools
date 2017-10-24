#include <vector>
#include <tuple>
#include <algorithm>
#include <math.h>
#include "perf.h"
#include "SeisppError.h"
#include "dmatrix.h"
#include "VectorStatistics.h"
#include "Vector3DBootstrapError.h"
using namespace SEISPP;
/* Important - through routine assumes x vectors are unit vectors.   Perhaps should verify this, but
   for efficiency probably won't do that. */
Vector3DBootstrapError::Vector3DBootstrapError(dmatrix& x,
    const double confidence, const int number_trials)
{
  const string base_error("Vector3DBootstrapError constructor:  ");
  /* Sanity check */
  if((confidence>1.0) || (confidence<=0.0))
  {
    throw SeisppError(base_error
        + "Illegal confidence interval requested - must be probability level (i.e greater than 0 and less than 1.0)");
  }
  cl=confidence;
  int nx=x.columns();
  int i,j,k;
  /* This is a large memory algorithm.  We load the trials into this matrix */
  dmatrix trials(3,number_trials);
  int col;
  double *swrk=new double[nx];
  //DEBUG
  //cout << "Computing trials matrix"<<endl;
  for(i=0;i<number_trials;++i)
  {
    dmatrix resampled(3,nx);
    for(j=0;j<nx;++j)
    {
      col=random_array_index(nx);
      for(k=0;k<3;++k) resampled(k,j)=x(k,col);
    }
    for(k=0;k<3;++k)
    {
      dcopy(nx,resampled.get_address(k,0),3,swrk,1);
      VectorStatistics<double> rstat(swrk,nx);
      trials(k,i)=rstat.mean();
    }
  }
  delete [] swrk;
  /*  Now from the suite of resampled trail values we need to
  compute the boostrap average and errors for angles expressed
  as dot products.  First compute the average of averages */
  vector<double> work;
  work.reserve(number_trials);
  /* This initialization might not be necessary*/
  for(i=0;i<number_trials;++i) work.push_back(0.0);
  double med[3];
  double nrmmed(0.0);
  //DEBUG
  //cout << "Computing mean vector"<<endl;
  for(k=0;k<3;++k)
  {
    dcopy(number_trials,trials.get_address(k,0),3,&(work[0]),1);
    VectorStatistics<double> vswrk(work);
    med[k]=vswrk.mean();
    nrmmed += (med[k]*med[k]);
  }
  /* Normalize the average vector */
  nrmmed=sqrt(nrmmed);
  for(k=0;k<3;++k) med[k] /= nrmmed;
  this->mean.clear();
  for(k=0;k<3;++k) this->mean.push_back(med[k]);
  /* Compute angles from average by dot product.  We reuse
     work - no need to clear as the above loop ends that way*/
  //DEBUG
  //cout << "Computing array of dot products"<<endl;
  for(i=0;i<number_trials;++i)
  {
    double dotprod,theta_trial,nrmtrials;
    nrmtrials=dnrm2(3,trials.get_address(0,i),1);
    for(k=0,dotprod=0.0;k<3;++k)
    {
      dotprod += trials(k,i)*med[k]/nrmtrials;
    }
    theta_trial=acos(dotprod);
    work[i]=theta_trial;
  }
  sort(work.begin(),work.end());
  /* This angle error is one sided - we estimate the probability
     the uncertainty in theta angles is less than the
     confidence value */
  int nconf=rint(((double)number_trials)*confidence);
  //DEBUG
  //cout << "Computed position for confidence="<<nconf<<" of "<<number_trials<<endl;
  if(nconf==number_trials) nconf=number_trials-1;   // Silently handled for roundoff errors
  if(nconf>=number_trials)
    throw SeisppError(base_error + "Coding problem - computed position of confidence outside array bound");
  aci=work[nconf];
}
pair<double,double> bootstrap_mv(const vector<double>& x, const double ci, const double ntrials)
{
  try{
    const string base_error("bootstrap_mv:  ");
    /* Sanity check */
    if((ci>1.0) || (ci<=0.0))
    {
      throw SeisppError(base_error
          + "Illegal confidence interval requested - must be probability level (i.e greater than 0 and less than 1.0)");
    }
    vector<double> trials;
    trials.reserve(ntrials);
    int nx=x.size();
    double *work=new double[nx];
    int i,j;
    for(j=0;j<ntrials;++j)
    {
      for(i=0;i<nx;++i)
      {
        int col;
        col=random_array_index(nx);
        work[i]=x[col];
      }
      VectorStatistics<double> vswrk(work,nx);
      trials.push_back(vswrk.mean());
    }
    delete [] work;
    sort(trials.begin(),trials.end());
    //DEBUG
    /*
    cout << "Sorted bootstrap mean values"<<endl;
    for(i=0;i<ntrials;++i) cout << trials[i]<<endl;
    */
    int ilow,ihigh,imed;
    double ltail,htail;
    ltail=(1.0-ci)/2.0;  //probability of smaller than lower cofidence limit
    htail=ltail+ci;  //probability of larger outer confidence limit
    ilow=rint( ( (double)ntrials)*ltail);
    ihigh=rint( ( (double)ntrials)*htail);
    //debug
    /*
    cout << "ltail, ilow, htail, ihigh"<<endl
      << ltail <<" "<<ilow<<" "<<htail<<" "<<ihigh<<endl;
    cout << trials[ihigh]<<" "<<trials[ilow]<<endl;
    */
    imed=ntrials/2;  // Assume ntrials is large enough to not worry abound odd/even
    pair<double,double> result;
    result.first=trials[imed];
    result.second = (trials[ihigh]-trials[ilow])/2.0; // return half range aka sigma in normal distribibution
    return result;
  }catch(...){throw;};
}
