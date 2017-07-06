#ifndef _PMTimeSeries_h_ 
#define _PMTimeSeries_h_
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "MWTransform.h"
#include "ParticleMotionEllipse.h"
#include "ParticleMotionError.h"
class PMTimeSeries : public BasicTimeSeries, public Metadata
{
    public:
        /*! Default constructor - initializes everything to zeros*/
        PMTimeSeries();
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
        /*! \brief Return the ellipse by sample number.

          This object encapsulates the concept of time-variable
        particle motion parameterized by an ellipse measured on
        a regular grid in time.   This method returns the ellipse
        by a time series type index.

        \param i sample number to return.
        \exception Throws a SeisppError object if i is out of range. 
        */
        ParticleMotionEllipse& ellipse(int i);
        /*! \brief Return the error statistics by sample number.

          This object encapsulates the concept of time-variable
        particle motion parameterized by an ellipse measured on
        a regular grid in time.   This method returns an object
        that contains all error estimates computed by the 
        multiwavelet tranform from a 3C seismogram.

        \param i sample number to return.
        \exception Throws a SeisppError object if i is out of range. 
        */
        ParticleMotionError errors(int i);


        vector<ParticleMotionEllipse> get_pmdata();
        vector<ParticleMotionError> get_pmerr();


        /*! Zero any data defined by a gap. 

          Warning this is not fully implemented. */
        void zero_gaps();
        /*! Standard assignment operator */
        PMTimeSeries& operator=(const PMTimeSeries& parent);
        /*! \brief Output data as ascii text.

          The intent of this procedure is to output the data from
          this object in a format that could be read to reconstruct
          it eventually.  At present it is a debug routine that
          just dumps the contents in a readable form.
          */
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
        friend class boost::serialization::access;
        template<class Archive>
                void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Metadata>(*this);
            ar & boost::serialization::base_object<BasicTimeSeries>(*this);
            ar & pmdata;
            ar & pmerr;
            ar & f0;
            ar & fw;
            ar & decfac;
        };
};
#endif
