#ifndef _PMTimeSeries_h_ 
#define _PMTimeSeries_h_
#include "ComplexTimeSeries.h"
#include "MWTransform.h"
#include "ParticleMotionEllipse.h"
#include "ParticleMotionError.h"
class PMTimeSeries : public BasicTimeSeries, public Metadata
{
    public:
        /*! Construct sample by sample for a specified band. 

          This will compute particle motions on a sample by 
          sample basis averaging overly over the bank of 
          multiwavelets.  This produces the maximum resolution
          and is the faster method to compute.  

          \param d is the input data computed with the multiwavelet
             transform and formed into a 3C bundle.
          \param band is the band from which this is to be computed.

          \exception SeisppError can be thrown for several illegal
             conditions. */
        PMTimeSeries(MWTBundle& d, int band);
        /*! Construct for a specified band with time average.

          This will compute particle motions with a dual averaging
          scheme.   First, the best fit ellipse is computed using 
          the largest singular vector of a matrix formed from each
          multiwavelet bundle from a common band over a specied 
          averaging length.   That is repeated for each wavelet
          and the particle motions are averaged to produce the 
          final estimate.   This process is repeated at a specified
          decimation interval over the length of the data.
          Note this is a vastly more expensive calculation than
          the sample by sample method with quadratic scaling 
          (theoretically) by the averaging length. 

          \param d is the input data computed with the multiwavelet
             transform and formed into a 3C bundle.
          \param band is the band from which this is to be computed.
          \param timestep is a decimation factor in initial sample
            interval steps.  e.g. if 4 the output sample interval
            with be 4*dt given dt as the input sample interval.
          \param avlen is the number of samples for each time average.

          \exception SeisppError can be thrown for several illegal
             conditions. */
        PMTimeSeries(MWTBundle& d, int band,int timesteps, int avlen);
        /*! Standard copy constructor. */
        PMTimeSeries(const PMTimeSeries& parent);
        ParticleMotionEllipse& ellipse(int i);
        ParticleMotionError errors(int i);
        void zero_gaps();
        PMTimeSeries& operator=(const PMTimeSeries& parent);
        friend ostream& operator<<(ostream& os, PMTimeSeries& d);
    private:
        vector<ParticleMotionEllipse> pmdata;
        vector<ParticleMotionError> pmerr;
        /* These attributes are stored here but should be 
           posted to Metadata to simplify interface. */
        int averaging_length;
        double f0,fw;
        int decfac;
        void post_attributes_to_metadata();
};
#endif
