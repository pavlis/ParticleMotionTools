#ifndef _ParticleMotionEllipse_h_
#define _ParticleMotionEllipse_h_
#include "TimeWindow.h"
#include "ComplexTimeSeries.h"
#include "dmatrix.h"
#include "MWTransform.h"
using namespace std;
using namespace SEISPP;
/*! \brief Particle motion ellipse computed by Multiwavelet method.

  The multiwavelet method allows the computation of particle motion ellipses
  from three-component seismic data.  The ellipse is defined by a rotational
  period defined by the center frequency of the analysis bandwidth.   This 
  object captures that concept.

  It has some oddities because of the need to interface it to older C code
  that I did not wish to mess with since it was know to work.   The big
  thing is that the attributes are all public to make it easier to deal with
  this interfacing. */
class ParticleMotionEllipse
{
public:
    /*! Unit vector defining the major axis direction in 3 space */
    double major[3];
    /*! Unit vector defining the minor axis direction in 3 space */
    double minor[3];
    /*! Size of major axis. */
    double majornrm;
    /*! Size of minor axis. */
    double minornrm;
    /*! Default constructor - initializes all attributes to zero. */
    ParticleMotionEllipse();
    /*! \brief Construct with a time averaging method.

      This method constructs a ParticleMotionEllipse from the output of
      the multiwavelet in a specified time window.   The method uses
      is a principal component method using a complex valued singular
      value decomposition.   

      \param x component in the x1 direction (normally +east)
      \param y component in the x2 direction (normally +north)
      \param z component in the x3 direction (normally +up)
      \param w defines the range in time to compute particle motion ellipse
         from the x,y,z data. 
      \param up is used to define the up direction.   This matters because there
       is a sign ambiguity in the sense of the vector that defines the
       maximum and minimum axes.  We resolve the ambituity by using the 
       direction where the dot product with the up vector is positive.  
         */

    ParticleMotionEllipse(ComplexTimeSeries& x, 
            ComplexTimeSeries& y, ComplexTimeSeries& z,
            TimeWindow w,double up[3]);
    /*! \brief Construct from complex numbers.

      A particle motion ellipse can be defined uniquely in 3-space from
    three complex numbers that define the amplitude and phase of a harmonic
    that makes exactly one cycle in 2*M_PI added to the phase.   This
    method uses an analytic formula to solve for the ellipse defined 
    by these three numbers. 

    \param x complex amplitude for component in x1 direction (normally +east)
    \param y complex amplitude for component in x2 direction (normally +north)
    \param z complex amplitude for component in x3 direction (normally +up)
    \param up is used to define the up direction.   This matters because there
       is a sign ambiguity in the sense of the vector that defines the
       maximum and minimum axes.  We resolve the ambituity by using the 
       direction where the dot product with the up vector is positive.  
       */
    ParticleMotionEllipse(Complex x,Complex y,Complex z,double up[3]);
    /*! \brief Construct from major and minor axis vectors.

      A particle motion ellipse is uniquely defined by its pricipal 
      component vector major an minor axes. This is the simplest 
      constructor that simply normalizes this pair of vectors and
      computes their lengths. */
    ParticleMotionEllipse(double *major,double *minor);
    /*! Standard copy constructor. */
    ParticleMotionEllipse(const ParticleMotionEllipse& parent);
    double rectilinearity();
    double major_azimuth();
    double major_inclination();
    double minor_azimuth();
    double minor_inclination();
    /*! Produce points to draw the ellipse in 3 space.

      In visualization we need to draw a particle motion ellipse
      in 3d space.   This method builds a matrix npointsx3 matrix of
      points to define the ellipse.

      \param npoints - number of points used to define the ellipse
      (Note there are no duplicate points.  To build a closed polygon repeat
      the first point after that last.)
      */
    dmatrix points(int npoints);
    /*! Initialize all attributes to zero. */
    void zero();

    /*! Standard assignment operator. */
    ParticleMotionEllipse& operator=(const ParticleMotionEllipse& parent);
    /*! Output stream operator.

      Produces a single ascii line with major(1 to 3) and minor (1 to 3) scaled by
      amplitudes (i.e. 6 columns of ascii data with blank separators. */
    friend ostream& operator<<(ostream& os, ParticleMotionEllipse& pme);
};

#endif
