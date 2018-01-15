#include "PMTimeSeries.h"
/* \brief Return an angle standardized to -180 to 180 degrees.
 *
 * In particle motion analysis there are lots of ways for 
 * computed angles to suffer wraparound problems.   We thus need ways
 * to standardize angles.   This functions standardizes any input angle (phi)
 * to be from -180 to 180 degrees.  Optional do same in radian units.
 *
  \param phi - input angle
  \param radians - if true assume phi is radians (default assumes degrees)
 
  \throw SeisppError if the input angle is absurd.
  \return angle in standard range
  */
double regularize_angle(double phi,bool radians)
{
    /* This is a necessary sanity check to avoid large unnecessary effort 
     * reducing an absurd number to standard range*/
    if( fabs(phi) > 100000.0) throw SeisppError(string("regularize_phi: ")
            + "Function was passed an absurd angle - coding error");
    if(radians)
    {
        if(phi>M_PI)
        {
            while(phi>M_PI)
            {
                phi -= 2.0*M_PI;
            }
        }
        else if(phi<(-M_PI))
        {
            while(phi<M_PI)
            {
                phi += 2.0*M_PI;
            }
        }
    }
    else
    {
        if(phi>180.0)
        {
            while(phi>180.0)
            {
                phi -= 360.0;
            }
        }
        else if(phi<(-180.0))
        {
            while(phi<180.0)
            {
                phi += 360.0;
            }
        }
    }
    return phi;
}
