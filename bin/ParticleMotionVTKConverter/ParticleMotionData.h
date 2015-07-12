/*! \brief This object encapsulates data for drawing 3C particle motion figures.

  This is the core data object for drawing three-component seismogram figures
  in a 3D visualization.  Concepts were worked out in PhD theses by 
  Anderson and Repin in the 1990s.   
  */
#include <vector>
#include "TimeWindow.h"
#include "ThreeComponentSeismogram.h"
#include "RegionalCoordinates.h"
#include "dmatrix.h"
class ParticleMotionData : public ThreeComponentSeismogram
{
    public:
        ParticleMotionData(ThreeComponentSeismogram& d, 
                RegionalCoordinates& x0);
        ParticleMotionData(const ParticleMotionData& parent);
        /*! Cartesian coordinates of center of figure to be drawn. */
        vector<double> center()
        {
            vector<double> result;
            result.reserve(3);
            for(int i=0;i<3;++i) result.push_back(x0[i]);
            return(result);
        };
        /*! Geographic coordinates of center of figure to be drawn. */
        Geographic_point geographic_center()
        {
            return x0_geo;
        };
        /*! \brief Return 3xN matrix of points defining the particle motion figure.
          \param tw is the time window to extract for this figure. 
          \exception SeisppException thrown if the window does not intersect
             the data range. */
        dmatrix particle_motion(TimeWindow tw);
        /*! \brief Return 3xN matrix of points defining the particle motion figure.
          This function returns points for a figure defined by the entire
          range of data stored in the object. */
        dmatrix particle_motion();
        /*! \brief Return 3xN matrix of raw data points.

          \param tw is the time window to extract for this figure. 
          \exception SeisppException thrown if the window does not intersect
             the data range. */
        dmatrix raw_motion(TimeWindow tw);
        dmatrix raw_motion();
        void set_scale(double sc)
        {
            scale=sc;
        };
        double get_scale(){return(scale);};
        /* These return maximum scaled amplitude */
        double max_amplitude();
        double max_amplitude(TimeWindow tw);
        /* Returns time window of full time span for data store here */
        TimeWindow timespan()
        {
            return(TimeWindow(this->t0,this->endtime()));
        };
        ParticleMotionData& operator=(const ParticleMotionData& parent);
    private:
        /* This determines scale factor to convert raw 3c data to positions
           in virtual output space.  Defaults to 1.0, which is 
        rarely appropriate. */
        double scale; 
        /* Data for the object are stored as geographic coordinates.
           This is the converter that converts to cartesian system.*/
        RegionalCoordinates coords;
        /* Coordinates of point defining the zero position for this 
           figure.  Cartesian coordinates that link to coords object.*/
        double x0[3];
        /* Geographic point equivalent to x0.   Cached for efficiency
           but use only for robustness.   This has to be extraceted from
           parent 3c seismogram metadata.   Idea is that setting it 
           internally means it must be defined and correct. */
        Geographic_point x0_geo;
        /* converts raw sample data to particle motion track - called by
           core method */
        void cook(dmatrix& d);
};


