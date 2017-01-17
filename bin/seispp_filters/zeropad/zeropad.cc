#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "seispp.h"
#include "ensemble.h"
using namespace std;   // most compilers do not require this
using namespace SEISPP;  //This is essential to use SEISPP library
ThreeComponentSeismogram pad_3cseis(ThreeComponentSeismogram& d,double plen,double tlen)
{
    int i,j;
    /* Return the original if marked dead */
    if(!d.live) return(d);
    /* Sanity checks on plen and tlen */
    const string base_error("pad_3cseis procedure:  ");
    double tracelength=d.endtime()-d.t0;
    if(plen>tracelength)
        throw SeisppError(base_error + "pad length irrational - exceeds data time length");
    const double taper_max_fraction(0.5);
    if(tlen>taper_max_fraction*tracelength)
        throw SeisppError(base_error + "taper length irrational - exceeds half of data time span");
    int nsnew,npad,ntaper;
    npad=plen/d.dt;
    ntaper=tlen/d.dt;
    nsnew=d.ns+npad;
    try {
        /* This buffer holds new data */
        dmatrix upad(3,nsnew);
        double *uptr;
        uptr=upad.get_address(0,0);
        /* This assumes dmatrix uses a single block of memory to hold matrix*/
        for(i=0;i<3*npad;++i,++uptr) *uptr=0.0;
        /* First copy the data before tapering */
        for(j=0;j<d.ns;++j)
            for(i=0;i<3;++i)
                upad(i,j+npad)=d.u(i,j);
        /* now do a linear taper */
        double wt0,wt;
        wt0=1.0/((double)ntaper);
        for(j=0,wt=wt0;j<ntaper-1;++j,wt+=wt0)
            for(i=0;i<3;++i)
                upad(i,j+npad)*=wt;
        d.u=upad;
        d.t0 -= plen;
        return d;
    }catch(...){throw;};
}

/* You shoudl always include a usage procedure like this to trap
   command line parsing problems. */
void usage()
{
    cerr << "zeropad [-pad dt -taper dt] < in > out"
        <<endl
        << "Zeropad all elements of a 3C ensemble"<<endl
        << " -pad sets zero pad time to dt (default 0.1)"<<endl
        << " -taper sets taper length to dt (default 0.01)"<<endl;
    exit(-1);
}
/* This is a generic routine to read a single object from 
   a boost text archive.   The object type is defined by
   the type InputObject as the template argument in standard
   C++ convention.   Note we return the object.   If 
   the object you are reading is huge you should consider 
   modifying this to return a pointer or auto_ptr.*/
template <class InputObject> InputObject 
    read_object(boost::archive::text_iarchive& ia)
{
    InputObject d;
    try{
        ia>>d;
    }catch(...)
    {
        /* This template ignores errors thrown by boost 
           and converts to a standard error for the seispp 
           library - perhaps not ideal but what is done here. */
        throw SeisppError(string("read_object failed:  ")
                    + "Check that input file is a boost text archive file");
    }
    return d;
}
/* This generic function is the inverse of read_object.  It writes
an object of type ObjectType to the boost archive stream.  
The write will fail with an exception if serialization is not 
defined for OutputObject. */
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
/* This obnoxious external variable is a necessary evil to deal with 
   error logging in the SEISPP library. Your code will probably not link 
   without it.*/ 
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    /* Common variables for a program common appear here, but 
       C/C++ now allow later declarations that usually make 
       cleaner code.   We need this counter for the arg list below*/
    int i;
    /* As the name implies set this to the number of required 
       args */
    const int narg_required(0);
    double padlength(0.1);
    double taperlength(0.01);

    for(i=narg_required+1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-pad")
        {
            ++i;
            if(i>=argc)usage();
            padlength=atof(argv[i]);
        }
        else if(sarg=="-taper")
        {
            ++i;
            if(i>=argc)usage();
            taperlength=atof(argv[i]);
        }
        else
            usage();
    }
    try{
        /* This allows input and output of objects through
           the boost serialization mechanism.   This assumes
           input and output of streams of data.  It has only
           been tested on single objects, but in principle it
           should work for arbitrarily complex formats mixing
           text data with boost archive text data. */
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentEnsemble d;
        d=read_object<ThreeComponentEnsemble>(ia);
        vector<ThreeComponentSeismogram>::iterator dptr;
        for(dptr=d.member.begin();dptr!=d.member.end();++dptr)
        {
            ThreeComponentSeismogram dpadded=pad_3cseis(*dptr,padlength,taperlength);
            *dptr=dpadded;
        }
        write_object<ThreeComponentEnsemble>(d,oa);
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(std::exception& stexc)
    {
        cerr << stexc.what()<<endl;
    }
}

