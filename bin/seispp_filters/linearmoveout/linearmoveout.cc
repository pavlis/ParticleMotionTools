#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "perf.h"
#include "seispp.h"
#include "ensemble.h"
/* This program applies a linear moveout computed from the 
   offset Metadata attribute stored with the data.  This 
   program is a unix filter that assumes input is a boost
   serialized ThreeComponentSeismogram object.   Output is 
   an altered version of the input.  

   This algorithm is less general than it should be because of
   a current limitation in the boost library.   I'm unable to
   get the boost library to save gap data.  Hence, here I 
   use the classic seismic processing approach of silently 
   zeroing gap data.  This will happen a lot with this algorithm.
   Whenever a negative t0 value is specified a zero pad will 
   always be created when the input is shot data.   If the input 
   is derived from continuous data that is a user error but it is 
   unavoidable with shot data.   

   This is my first attempt at a seispp processing unix style
   filter with boost serialization.   If this works I expect to 
   use the skeleton of this code to construct a template for 
   a "module" meshing with that concept.

Author:  Gary L Pavlis
Written:  Nov 29, 2015
*/
void usage()
{
    cerr << "linearmoveout [-vr v_reduce -t0 timeshift --help] < infile > outfile"
        <<endl
        << " -vr sets reducing velocity (default 6000)"<<endl
        << " -t0 applies a time shift to all seismograms (default 0)"<<endl;;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    double vreduce(6000.0);   // unit so m/s for this default
    double t0(0.0);
    int i;
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-t0")
        {
            ++i;
            if(i>=argc)usage();
            t0=atof(argv[i]);
        }
        else if(sarg=="-vr")
        {
            ++i;
            if(i>=argc)usage();
            vreduce=atof(argv[i]);
        }
        else if(sarg=="--help")
            usage();
        else
            usage();
    }
    try{
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentEnsemble d;
        ia >> d;
        int nm=d.member.size();
        for(i=0;i<nm;++i)
        {
            if(d.member[i].live)
            {
               ThreeComponentSeismogram work(d.member[i]);
               /* Compute the time shift for this signal as
                  t0+offset/vreduce  */
               double tshift=t0;
               double offset=work.get_double("offset");
               offset=fabs(offset);
               tshift -= (offset/vreduce);
               int ioffset=work.sample_number(tshift);
               /* a basic sanity check on velocity and t0  */
               if(abs(ioffset)>work.ns)
               {
                   cerr << "linearmoveout(Fatal Error):  irrational time shift"
                        <<endl;
                   cerr << "using t0="<<t0<<" and reducing velocity="<<vreduce
                       << endl
                       << "This seismogram has offset set as "<<offset
                       << " which yields a time shift of "<<ioffset<<"samples"
                       <<endl
                       << "This exceeds number of data samples ="
                       <<work.ns<<endl;
                   exit(-1);
               }
               int nsout=work.ns - ioffset;
               int ncopy=min(work.ns,nsout);
               work.u=dmatrix(3,nsout);
               work.u.zero();
               work.t0+=t0;
               int lag=work.sample_number(-offset/vreduce);
               //DEBUG
               cerr << "lag="<<lag
                   << " nsout="<< nsout 
                   << " offset="<< offset 
                   << " ioffset="<< ioffset 
                   << " ncopy="<< ncopy 
                   <<endl;
               double *ptr;
               ptr=work.u.get_address(0,lag);
               dcopy(3*ncopy,d.member[i].u.get_address(0,0),1,ptr,1);
               work.ns=nsout;
               d.member[i]=work;
            }
        }
        oa << d;
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(std::exception& stexc)
    {
        cerr << stexc.what()<<endl;
    }
}

