#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "perf.h"
#include "seispp.h"
#include "ensemble.h"
using namespace std;   // most compilers do not require this
using namespace SEISPP;  //This is essential to use SEISPP library
void usage()
{
    cerr << "sphdiv [-decay power] < infile >outfile"
        <<endl
        << "Applies spherical divergence correction to three component"<<endl
        << "data stored in a boost test archive file."<<endl
        << "Amplitudes scaled by metadata offset variable to power."<<endl
        << "There is no normalization so beware"<<endl
        << "Default power is 2.0"<<endl;
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
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    int i;
    const int narg_required(0);
    double decay_power(2.0);

    for(i=narg_required+1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-power")
        {
            ++i;
            if(i>=argc)usage();
            decay_power=atof(argv[i]);
        }
        else
            usage();
    }
    try{
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentEnsemble d;
        d=read_object<ThreeComponentEnsemble>(ia);
        vector<ThreeComponentSeismogram>::iterator dptr;
        double scale;
        const double MIN_offset(0.001);
        for(dptr=d.member.begin();dptr!=d.member.end();++dptr)
        {
            /* Do nothing if the seismogram is marked bad */
            if(dptr->live)
            {
                double offset=dptr->get_double("offset");
                /* sanity check to avoid zeroing data.  
                   Assumes rational units. */
                if(offset>MIN_offset)
                {
                    scale=pow(offset,decay_power);
                    /* ThreeComponentSeismogram really needs an operator
                       that does this kind of scale - operator * */
                    dscal(3*dptr->ns,scale,dptr->u.get_address(0,0),1);
                }
            }
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

