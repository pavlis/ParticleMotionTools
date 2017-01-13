#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "seispp.h"
#include "ensemble.h"
using namespace std;
using namespace SEISPP;

void usage()
{
    cerr << "extract_component n < infile > outfile"<<endl
        << "Extract component n (must be 0, 1, or 2) from input 3C ensemble"
        <<endl
        << "Output is boost serialized TimeSeriesEnsemble object"<<endl;
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

bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{

    int i;
    const int narg_required(1);
    if(argc!=narg_required) usage();
    int outchan=atoi(argv[1]);
    if( (outchan<0) || (outchan>2) )
    {
      cerr << "extract_component:  illegal output channel of "
        <<outchan<<" specified"<<endl;
      usage();
    }
    try{
        if(SEISPP_verbose)
        {
          cerr << "extract_component:  extracting component number "
            << outchan<<endl;
        }
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentEnsemble d;
        d=read_object<ThreeComponentEnsemble>(ia);
        auto_ptr<TimeSeriesEnsemble> dscalar=ExtractComponent(d,outchan);
        write_object<TimeSeriesEnsemble>(*dscalar,oa);
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(std::exception& stexc)
    {
        cerr << stexc.what()<<endl;
    }
}
