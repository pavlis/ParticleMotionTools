/* This is a filter that will separate a TimeSeries of ThreeComponent
   seismogram into components stripping ensemble header values and
   putting copies into the header of each member.   This is a required
   step before sort or any other process that operates on one seismogram
   at a time.
*/

/* This set of system includes are always required.  Do not remove them.*/
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "seispp.h"
#include "ensemble.h"
using namespace std;   // most compilers do not require this
using namespace SEISPP;  //This is essential to use SEISPP library
void usage()
{
    cerr << "dismember [-scalar]"
        <<endl
        << "seispp filter separates ensemble into seismogram members"<<endl
        << " Use -scalar for TimeSeriesEnsembles (default is ThreeComponentEnsemble data)"
        <<endl;
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
/* Generic algorithm to serialize an ensemble.   Copies the
ensemble metadatra to each output seismogram.   Returns
count of the number of seismograms written.*/
template <class Tens,class Tmem> int write_ensemble(Tens& d,
        boost::archive::text_oarchive& oa)
{
    try {
      int i;
      Metadata ensmd(dynamic_cast<Metadata&>(d));
      MetadataList keylist=ensmd.keys();
      for(i=0;i<d.member.size();++i)
      {
        /* this little procedure makes the copy of metadata
        pretty trivial */
        copy_selected_metadata(ensmd,dynamic_cast<Metadata&>(d.member[i]),keylist);
        oa << d.member[i];
      }
      return i;
    }catch(...)
    {
        throw SeisppError(string("write_object failed\n")
                +"Is serialization defined for this object type?\n"
                +"Do you have write permission for output directory?");
    }
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    int i;
    bool scalar_mode(false);
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-scalar")
          scalar_mode=true;
        else
            usage();
    }
    try{
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentEnsemble d3c;
        TimeSeriesEnsemble d;
        int nensembles(0),nd(0);
        /* This approach is the one described at this URL to
           stop on EOF - not the most elegant since it requires
           one to use an exception to terminate the output, but
           if it works I'll accept that.
http://stackoverflow.com/questions/7111041/boost-serialization-multiple-objects
        */
        try {
            for(;;)
            {
                int count;
                if(scalar_mode)
                {
                  d=read_object<TimeSeriesEnsemble>(ia);
                  count=write_ensemble<TimeSeriesEnsemble,TimeSeries>
                      (d,oa);
                  nd+=count;
                }
                else
                {
                  d3c=read_object<ThreeComponentEnsemble>(ia);
                  count=write_ensemble<ThreeComponentEnsemble,ThreeComponentSeismogram>
                      (d3c,oa);
                  nd+=count;
                }
                ++nensembles;
            }
        }catch(boost::archive::archive_exception const& e)
        {
          cerr << "Read "<<nensembles<<" and wrote "<<nd<<" seismograms to output"<<endl;
          cerr << "boost archive error message used to catch eof"<<endl;
          cerr << e.what()<<endl;
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
