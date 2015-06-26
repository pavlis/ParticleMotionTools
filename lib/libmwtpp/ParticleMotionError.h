#include <iostream>
using namespace std;
/* \brief Data object to hold multiwavelet generated error estimates.

   This object holds attributes that are computed from multiwavelet 
   analyis of three-component seismogram data to yield a particle 
   motion ellipse estimate.  Multiwavelets have redundant estimates
   that can be used to form these error estimates.

   This object is essentially a C struct converted to a C++ object.
   It is largely used to hold output of an older C procedure.*/
class ParticleMotionError
{
public:
    /*! Error in emergence angle of major axis of ellipse. */
    double dtheta_major; 
    /*! Error in azimuth angle of major axis of ellipse. */
    double dphi_major;
    /*! Error in emergence angle of minor axis of ellipse. */
    double dtheta_minor; 
    /*! Error in azimuth angle of minor axis of ellipse. */
    double dphi_minor;
    /*! Error in size of major axis of ellipse. */
    double dmajornrm; 
    /*! Error in size of minor axis of ellipse. */
    double dminornrm;
    /*! Error in rectilinearity parameter of ellipse. */
    double delta_rect;
    /*! Number of degrees of freedom of major axis parameters. */
    int ndgf_major; 
    /*! Number of degrees of freedom of minor axis parameters. */
    int ndgf_minor; 
    /*! Number of degrees of freedom of rectilinearity estimate. */
    int ndgf_rect;
    /*! Number of degrees of freedom of major axis amplitude estimate. */
    int ndgf_major_amp; 
    /*! Number of degrees of freedom of minor axis amplitude estimate. */
    int ndgf_minor_amp;
    /*! Defaault constructor.   

      Initializes all data to zero. */
    ParticleMotionError();
    /*! Zero all attributes - convenience function. */
    void zero();
    /*! Output stream operator.

      The output is assumed to be in ascii.  Output is one line with
      no newline character at the end.   The order is the same as the 
      list of attributes in the class above. */
    friend ostream& operator<<(ostream& os, ParticleMotionError& pme);
};
