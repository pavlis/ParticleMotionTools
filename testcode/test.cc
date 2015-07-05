#include <stdlib.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "seispp.h"
#include "MWTransform.h"
#include "PMTimeSeries.h"
using namespace std;
using namespace SEISPP;
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    try {
        Complex xs(1.0,0.0);
        Complex ys(2.0,0.0);
        Complex zs(0.0,0.0);
        double up[3]={0.0,0.0,1.0};
        ParticleMotionEllipse pm1(xs,ys,zs,up);
        cout << "Linear polarized test"<<endl;
        cout << "Rectilinearity="<<pm1.rectilinearity()<<endl;
        cout << "PM ellipse contents:"<<endl;
        cout <<pm1<<endl;
        cout << "Points"<<endl;
        dmatrix p=pm1.points(360);
        cout << p<<endl;
        cout << "Trying to use serialize )code for PMEllipse"<<endl<<"Writing to ellipse.dat"<<endl;
        std::ofstream ofs("ellipse.dat");
        boost::archive::text_oarchive oa(ofs);
        oa << pm1;
        ofs.close();
        cout << "Success:  Trying to restore."<<endl;
        std::ifstream ifspe("ellipse.dat");
        boost::archive::text_iarchive ia(ifspe);
        ParticleMotionEllipse pmread;
        ia >> pmread;
        ifspe.close();
        cout << "Read completed.  Contents read back follow:"<<endl;
        cout << pmread;
        ThreeComponentSeismogram d(2000);
        d.ns=2000;
        d.t0=0.0;
        d.dt=0.01;
        d.tref=relative;
        d.live=true;
        d.put("samprate",100.0);
        d.put("time",0.0);
        d.put("nsamp",2000);
        d.u.zero();
	/* Initialize to small random numbers to avoid nans */
	int i,k;
	int offset=RAND_MAX/2;
	for(i=0;i<2000;++i)
		for(k=0;k<3;++k)
		{
			d.u(k,i)=(double) ((random()-offset)/((double)RAND_MAX));
		}
        d.u(0,400)=1000.0;
        d.u(1,400)=2000.0;
        d.u(0,800)=1000.0;
        d.u(1,820)=2000.0;
        d.u(0,1200)=1000.0;
        d.u(1,1220)=1000.0;
/*
	dmatrix uraw=tr(d.u);
	cout << "Input 3C data matrix"<<endl;
	cout <<uraw<<endl;
*/
        string fname("test.pf");
        PfStyleMetadata md=pfread(fname);
        cout << dynamic_cast<Metadata&> (md)<<endl;
	cout << "Creating multiwavelet transformer object"<<endl;
        MWTransform mwt(fname);
	cout << "Success"<<endl<<"Computing transform of test waveform"<<endl;
        MWTBundle dtrans(d,mwt);
	cout << "Success"<<endl<<"Data from band 0, wavelet 0, component 0"<<endl;
	MWTwaveform d00=dtrans(0,0,0);
	cout << dynamic_cast<ComplexTimeSeries&>(d00)<<endl;
        PMTimeSeries pmts(dtrans,0);
        cout << "Testing boost serialization output to file PMtest.dat"<<endl;
        std::ofstream ofs2("PMTtest.dat");
        boost::archive::text_oarchive opmts(ofs2);
        opmts << pmts;
        cout << "Save with serialization test write successed"<<endl
            << "Attempting inverse operation"<<endl;
        ofs2.close();
        std::ifstream ifs("PMTtest.dat");
        boost::archive::text_iarchive ipmts(ifs);
        PMTimeSeries pmtsin;
        ipmts >> pmtsin;
        cout << "Read completed:  "<<endl
            << "ns of data written="<<pmts.ns<<" ns of data read from file="<<pmtsin.ns<<endl
            <<"Comparing contents by testing major axis values only"<<endl
            << "major written (3 columns);  major read (3 columns);  difference (3 columns)"<<endl;
        for(i=0;i<pmts.ns;++i)
        {
            ParticleMotionEllipse pmo(pmts.ellipse(i));
            ParticleMotionEllipse pmi(pmtsin.ellipse(i));
            for(k=0;k<3;++k) cout << pmo.major[k]<<" ";
            for(k=0;k<3;++k) cout << pmi.major[k]<<" ";
            for(k=0;k<3;++k) cout << pmo.major[k]-pmi.major[k]<<" ";
            cout <<endl;
        }
        
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(std::exception& err)
    {
        cout << err.what()<<endl;
    }
}
