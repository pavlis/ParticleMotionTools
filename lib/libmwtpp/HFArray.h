#include <vector>
#include <map>
#include "slowness.h"
#include "RegionalCoordinates.h"
#include "dbpp.h"
using namespace SEISPP;
using namespace std;
class HFArray 
{
    public:
        /*! \brief Construct from a text file.

          */
        HFArray(string fname);
        /*! \brief Construct from a css3.0 database. 

          The antelope css3.0 schema database uses a site table to 
          define locations.   This constructor uses only the dnorth
          and deast attributes to define array coordinates. 
          
          \param dbh is a handle the Antelope db (does not need to 
             point at site
          \param net is SEED netname to limit db site selection.  If
          this is empty the entire site table will be loaded.
        */
        HFArray(DatascopeHandle& dbh,string net);
        /*! Standard copy constructor. */
        HFArray(const HFArray& parent);
        /*! Return Cartesian coordinates (in km) as e,n,z vector. */
        vector<double> x(string sta);
        /* Compute moveout for surface array */
        double moveout(string sta, SlownessVector u);
        /* Compute moveout for 3D array with large elevation variations */
        double moveout(string sta, SlownessVector u, double v0);
        Geographic_point origin();
        Geographic_point geographic_location(string sta);
        HFArray& operator=(const HFArray& parent);
    private:
        map<string,vector<double> > stations;
        RegionalCoordinates coords;
};
