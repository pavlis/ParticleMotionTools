#include <vector>
#include <map>
#include "slowness.h"
#include "RegionalCoordinates.h"
using namespace SEISPP;
using namespace std;
class HFArray 
{
    public:
        HFArray(string fname);
        HFArray(const HFArray& parent);
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
