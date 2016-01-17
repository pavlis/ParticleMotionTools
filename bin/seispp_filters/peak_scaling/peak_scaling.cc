/* This small program scales each 3c seismogram by peak amplitude.  */
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
/* These are SEISPP includes that are a base requirement 
   to use this template file */
#include "seispp.h"
/* Replace these as needed.  This one is here only to make
   this do nothing template compile for testing configuration */
#include "ensemble.h"
/* You will get lots of errors without these namespace
   declaration*/
using namespace std;   // most compilers do not require this
using namespace SEISPP;  //This is essential to use SEISPP library
/* You shoudl always include a usage procedure like this to trap
   command line parsing problems. */
void usage()
{
    cerr << "peak_amplitude < infile > outfile"
        <<endl
        << "Reads serialized 3c ensemble file and scales each "
        << "seismogram by peak 3C amplitude"<<endl;
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
    if(argc!=1) usage();
    try{
        /* This allows input and output of objects through
           the boost serialization mechanism.   This assumes
           input and output of streams of data.  It has only
           been tested on single objects, but in principle it
           should work for arbitrarily complex formats mixing
           text data with boost archive text data. */
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        /* This uses the generic template above to read a single
           object.  This example uses a ThreeComponentEnsemble 
           but any object with serialization define should work
           with this template. Note use of an auto_ptr for 
           efficiency with large objects*/
        ThreeComponentEnsemble d;
        d=read_object<ThreeComponentEnsemble>(ia);
        vector<ThreeComponentSeismogram>::iterator dptr;
        int i;
        for(dptr=d.member.begin(),i=0;dptr!=d.member.end();++dptr,++i)
        {
            double amp;
            ThreeComponentSeismogram *ptr;
            ptr=&(d.member[i]);
            amp=PeakAmplitude(ptr);
            double gain(1.0);
            try {
                gain=dptr->get_double("gain");
            }catch(MetadataGetError &mde)
            {
                if(SEISPP_verbose)
                    cerr << "Warning:   gain attribute was not set, using default 1"
                        <<endl;
            }
            double scaling=1.0/amp;
            gain = gain*scaling;
            dptr->put("gain",gain);
            dptr->u=scaling*dptr->u;
            if(SEISPP_verbose)
                cerr << "Computed gain scaling="<<scaling<<endl;
        }
        /*This is the output equivalent of read_object.  This 
         is again only an example and should be changed.*/
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

