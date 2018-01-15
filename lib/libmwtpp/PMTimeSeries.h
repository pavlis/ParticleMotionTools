#ifndef _PMTimeSeries_h_ 
#define _PMTimeSeries_h_
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "MWTransform.h"
#include "BasicTimeSeries.h"
#include "Metadata.h"
#include "ParticleMotionEllipse.h"
#include "ParticleMotionError.h"
using namespace SEISPP;
/*! \brief Time series style representation of particle motion ellipse data.
 *
 The multiwavelet transform can be used to produce particle motion estimates
 as a function of time.   This object encapsulates that concept as an 
 abstraction of a TimeSeries object wherein the data points at each 
 time step are themselves objects = ParticleMotionEllipse implementation.
 */
/*! Metadata key that tags TimeSeries scalar metrics returned by
 * several methods in this object */
const string PMDerivedTSType("DerivedTSDataType");
/* Azimuth error estimates for particle motion data scale by
 * 1/sin(inclination).   When inclinations are estimated very small 
 * inclination estimates are common.   These two variables control
 * how this is handled in construction of PMTimeSeries objects.   
 * thetafloor is an absolute floor.   When the estimated inclination
 * angle is less than this number (below 1 degree in radians) the
 * azimuth errors are set to +-180 degrees (pi radians).  The second
 * floor is the error_thetafloor variable.  When the inclination angle
 * lies between thetafloor and error_thetafloor the azimuth error 
 * estimate is scaled by 1/sin(error_thetafloor).   We could code that
 * constant, but it would make this obscure for handling a rare event. 
 * i.e. the cost of computing 1/sin(floor) - i.e. cosecant - is small.
 * */
const double thetafloor(0.017453292519943);
const double error_inclination_floor(0.17453292519943);
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
          \param confidence specifies confidence level to compute using
             bootstrap (default 95%).
          \param bsmultiplier defines the number of trials to use in
             computing bootstrap errors 
             (number of trails=bsmultiplier*number_of_wavelets). 

          \exception SeisppError can be thrown for several illegal
             conditions. */
        PMTimeSeries(MWTBundle& d, int band, 
            double confidence=0.95,int bsmultiplier=100);
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
          \param confidence specifies confidence level to compute using
             bootstrap (default 95%).
          \param bsmultiplier defines the number of trials to use in
             computing bootstrap errors 
             (number of trails=bsmultiplier*number_of_wavelets). 

          \exception SeisppError can be thrown for several illegal
             conditions. */
        PMTimeSeries(MWTBundle& d, int band,int timesteps, int avlen,
            double confidence=0.95,int bsmultiplier=100);
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


        /*! Get the raw vector of particle motion data defined as ellipses.*/
        vector<ParticleMotionEllipse> get_pmdata();
        /*! Ge the raw vector of error estimates for the ellipses. */
        vector<ParticleMotionError> get_pmerr();
        /*! \brief Get representation of major axis length as a time series.

          A PMTimeSeries is used to save a full representation of 
          particle motion estimates as a function of time.  Various
          scalar metrics can be computed from particle motions 
          represented as an ellipse.   
          This method retrieves the L2 norm length of the major axis of 
          the particle motion ellipse at each time step.
          The Metadata attributes of the original PMTimeSeries are
          copied to the TimeSeries that is returned.  
          */
        TimeSeries major_axis_amplitude();
        /*! \brief Get representation of minor axis length as a time series.

          A PMTimeSeries is used to save a full representation of 
          particle motion estimates as a function of time.  Various
          scalar metrics can be computed from particle motions 
          represented as an ellipse.   
          the particle motion ellipse at each time step.
          The Metadata attributes of the original PMTimeSeries are
          The Metadata attributes of the original PMTimeSeries are
          copied to the TimeSeries that is returned.  
          */
        TimeSeries minor_axis_amplitude();
        /*! \brief Get representation of rectilinearity as a time series.

          A PMTimeSeries is used to save a full representation of 
          particle motion estimates as a function of time.  Various
          scalar metrics can be computed from particle motions 
          represented as an ellipse.   
          This method retrieves the metric called rectilinearity.   
          The Metadata attributes of the original PMTimeSeries are
          copied to the TimeSeries that is returned.  
          */
        TimeSeries rectilinearity();
        /*! \brief Get representation of major azimuth as a time series.

          A PMTimeSeries is used to save a full representation of 
          particle motion estimates as a function of time.  Various
          scalar metrics can be computed from particle motions 
          represented as an ellipse.   
          This method retrieves the major axis direction as function of
          time (in radians). 
          The Metadata attributes of the original PMTimeSeries are
          copied to the TimeSeries that is returned.  
          */
        TimeSeries major_azimuth();
        /*! \brief Get representation of major axis inclination  as a time series.

          A PMTimeSeries is used to save a full representation of 
          particle motion estimates as a function of time.  Various
          scalar metrics can be computed from particle motions 
          represented as an ellipse.   
          This method retrieves the major axis direction as function of
          time (in radians). 
          The Metadata attributes of the original PMTimeSeries are
          copied to the TimeSeries that is returned.  
          */
        TimeSeries major_inclination();
        /*! \brief Get representation of minor azimuth as a time series.

          A PMTimeSeries is used to save a full representation of 
          particle motion estimates as a function of time.  Various
          scalar metrics can be computed from particle motions 
          represented as an ellipse.   
          This method retrieves the minor axis direction as function of
          time (in radians). 
          The Metadata attributes of the original PMTimeSeries are
          copied to the TimeSeries that is returned.  
          */
        TimeSeries minor_azimuth();
        /*! \brief Get representation of minor axis inclination  as a time series.

          A PMTimeSeries is used to save a full representation of 
          particle motion estimates as a function of time.  Various
          scalar metrics can be computed from particle motions 
          represented as an ellipse.   
          This method retrieves the minor axis direction as function of
          time (in radians). 
          The Metadata attributes of the original PMTimeSeries are
          copied to the TimeSeries that is returned.  
          */
        TimeSeries minor_inclination();


        /*! Zero any data defined by a gap. 

          Warning this is not fully implemented. */
        void zero_gaps();
        /*! Standard assignment operator */
        PMTimeSeries& operator=(const PMTimeSeries& parent);
        double get_wavelet_duration(){return wavelet_duration;};
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
        /* Because we allow step the duration of the wavelet is stored
         * as a time interval. In MWTMatrix this is save by number of
         * samples */
        double wavelet_duration;
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
            /* This was an add on with significant amounts of data
             * stored without it.  This uses the method described
             * in the serialization tutorial to handle data versions.
             */
            if(version>0)
                ar & wavelet_duration;
        };
};
BOOST_CLASS_VERSION(PMTimeSeries,1);
double regularize_angle(double d,bool r=false);
#endif
