#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "seispp.h"
#include "ensemble.h"
#include "interpolator1d.h"
using namespace std;
using namespace SEISPP;

void usage()
{
    cerr << "agc windlength [-g gainfunction] < infile > outfile"
        <<endl
        << "Applies three-component AGC operator of duraction winlength second"
        <<endl
        << "Optionally save ensemble of gain function as TimeSeriesEnsemble  object with winlength/2 sample interval"
        <<endl;
    exit(-1);
}

template <class InputObject> InputObject
    read_object(boost::archive::text_iarchive& ia)
{
    InputObject d;
    try{
        ia>>d;
    }catch(...)
    {
        throw SeisppError(string("read_object failed:  ")
                    + "Check that input file is a boost text archive file");
    }
    return d;
}

template <class OutputObject> void write_object(OutputObject& d,
        boost::archive::text_oarchive& oa)
{
    try {
        oa << d;
    }catch(...)
    {
        throw SeisppError(string("write_object failed\n")
                +"Is serialization defined for this object type?\n"
                +"Do you have write permission for output directory?");
    }
}
TimeSeries agc(ThreeComponentSeismogram& d, double twin)
{
    try{
        TimeSeries gf(dynamic_cast<Metadata&>(d),false);
        gf.dt=twin/2.0;
        /* This puts the first sample of the gain function center of 
         * agc window */
        gf.t0=d.t0+gf.dt;
        gf.s.clear();
        /* Compute number of samples in agc window */
        int nsw=static_cast<int>((gf.dt)/(d.dt));
        nsw*=2;  // Assures nsw is devisible by 2 - require for loop below 
        if(nsw<=0) throw SeisppError("agc funtion:  window length is less than sample interval");
        double ssq,scale;
        int i,k,iw;
        for(iw=0;iw<d.ns;iw+=gf.dt)
        {
            for(i=iw;i<iw+nsw;++i)
                for(k=0;k<3;++k)
                {
                    double val=d.u(i,k);
                    ssq += val*val;
                }
            scale=sqrt(ssq)/((double)(3*nsw));
            gf.s.push_back(scale);
        }
        /* Now work through the data sample by sample and apply
         * the gain function.  We use a linear interpolator between
         * computed gain values */
        for(i=0;i<d.ns;++i)
        {
            //Reuse iw as position in gain function
            iw=gf.sample_number(d.time(i));
            /* This assumes sample_number uses nint so the point is between
             * iw and iw+1*/
            if(iw<0)
                scale=gf.s[0];
            else if(iw>=(d.ns-1))
                scale=gf.s[gf.ns-1];
            else
                scale=INTERPOLATOR1D::linear_scalar(d.time(i),
                        gf.time(iw),gf.s[iw],gf.time(iw+1),gf.s[iw+1]);
            for(k=0;k<3;++k) d.u(i,k)*=scale;
        }
        return gf;
    }catch(...){throw;};
}

bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{

    int i;
    const int narg_required(1);
    if(argc<2) usage();
    double agcwinlen=atof(argv[1]);
    bool save_gain_function(false);
    string gainfile;

    for(i=narg_required+1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-g")
        {
            ++i;
            if(i>=argc)usage();
            gainfile=string(argv[i]);
            save_gain_function=true;
        }
        else
            usage();
    }
    try{
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentEnsemble d;
        d=read_object<ThreeComponentEnsemble>(ia);
        int n=d.member.size();
        TimeSeriesEnsemble gains(dynamic_cast<Metadata&>(d),n);
        int i;
        for(i=0;i<n;++i)
        {
            TimeSeries g=agc(d.member[i],agcwinlen);
            gains.member.push_back(g);
        }
        write_object<ThreeComponentEnsemble>(d,oa);
        if(save_gain_function)
        {
            ofstream ofs;
            ofs.open(gainfile.c_str(),ios::out);
            if(ofs.fail())
            {
                cerr << "agc: Could not open file "<<gainfile
                    <<" to save gain functions for input ensemble"<<endl
                    << "Fix this if you need the gain function."<<endl
                    << "agc data written successfully to stdout"<<endl;
            }
            boost::archive::text_oarchive oagf(ofs);
            write_object<TimeSeriesEnsemble>(gains,oagf);
            ofs.close();
        }
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(std::exception& stexc)
    {
        cerr << stexc.what()<<endl;
    }
}
