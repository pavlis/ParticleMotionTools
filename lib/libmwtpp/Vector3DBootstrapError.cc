#include <vector>
#include <tuple>
#include <algorithm>
#include <math.h>
#include "SeisppError.h"
#include "dmatrix.h"
#include "VectorStatistics.h"
#include "Vector3DBootstrapError.h"
using namespace SEISPP;
/* Important - through routine assumes x vectors are unit vectors.   Perhaps should verify this, but
   for efficiency probably won't do that. */
Vector3DBootstrapError::Vector3DBootstrapError(dmatrix x,double confidence, int number_trials)
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
  for(i=0;i<number_trials;++i)
  {
    int col;
    col=random_array_index(nx);
    for(k=0;k<3;++i) 
       trials(k,i)=x(k,col);
  }
  vector<double> work;
  work.reserve(number_trials);
  med.reserve(3);
  double nrmmed(0.0);
  for(k=0;k<3;++k)
  {
    for(i=0;i<number_trials;++i)
    {
      work.push_back(trials(k,i));
    }
    med.push_back(median<double>(work));
    nrmmed += (med[k]*med[k]);
    work.clear();
  }
  /* Normalize the median vector */
  nrmmed=sqrt(nrmmed);
  for(k=0;k<3;++k) med[k] /= nrmmed;
  /* Compute angles from median by dot product.  We reuse
     work - no need to clear as the above loop ends that way*/
  for(i=0;i<number_trials;++i)
  {
    double dotprod,theta;
    for(k=0,dotprod=0.0;k<3;++k)
    {
      dotprod += trials(k,i)*med[k];
    }
    theta=acos(dotprod);
    work.push_back(theta);
  }
  sort(work.begin(),work.end());
  /* This angle error is one sided - we estimate the probability
     the uncertainty in theta angles is less than the
     confidence value */
  int nconf=rint(((double)number_trials)*confidence);
  if(nconf==number_trials) nconf=number_trials-1;   // Silently handled for roundoff errors
  if(nconf>=number_trials)
    throw SeisppError(base_error + "Coding problem - computed position of confidence outside array bound");
  aci=work[nconf];
}
