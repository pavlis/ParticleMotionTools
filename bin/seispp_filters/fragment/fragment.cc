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
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "stock.h" // needed for makedir call
#include "seispp.h"
#include "ensemble.h"
using namespace std;   // most compilers do not require this
using namespace SEISPP;  //This is essential to use SEISPP library
void usage()
{
    cerr << "fragment basename [-dir outdir -v]"
        <<endl
        << "seispp filter fragments file with multiple ensembles into individual files"
        <<endl
        << "basename is root name for each ensemble.  Adds a sequence number for each ensemble"
        <<endl
        << "-dir optional write to outdir (default is .)"
        <<endl
        << "-v verbose output (mostly logs each ensembles gather metadata"
        <<endl;
    exit(-1);
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
      }
      oa << d;
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
    bool Verbose(false);
    int i;
    if(argc<2) usage();
    string outdir(".");
    string basename(argv[1]);
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-dir")
        {
            ++i;
            if(i>=argc)usage();
            outdir=string(argv[i]);
            if(makedir(const_cast<char *>(outdir.c_str())))
            {
                cerr << "Cannot create requested directory "
                    << outdir<<endl;
                usage();
            }
        }
        else if(sarg=="-v")
            Verbose=true;
        else
            usage();
    }
    try{
        boost::archive::text_iarchive ia(cin);
        ThreeComponentEnsemble d3c;
        int nensembles(0);
        int nseis(0);
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
                ia >> d3c;
                char fname[128];
                sprintf(fname,"%s_%d",basename.c_str(),nensembles);
                string path;
                path=outdir+"/"+fname;
                if(Verbose)
                {
                    cerr << "ensemble (gather) metadata for output file "
                        << path<<endl;
                    cerr << dynamic_cast<Metadata&>(d3c)<<endl;
                }
                ofstream out(path.c_str(),ios::out);
                if(!out)
                {
                    cerr << "Cannot open output file " << path
                        <<endl;
                    usage();
                }
                boost::archive::text_oarchive oa(out);
                count=write_ensemble<ThreeComponentEnsemble,ThreeComponentSeismogram>
                    (d3c,oa);
                out.close();
                ++nensembles;
                nseis+=count;
            }
        }catch(boost::archive::archive_exception const& e)
        {
          if(Verbose){
          cerr << "Processed "<<nensembles<<" ensembles"<<endl
              << "with a total of "<<nseis<<" three component seismograms"<<endl
              << " Look for output in directory "<<outdir<<endl;
          cerr << "boost archive error message used to catch eof"<<endl;
          cerr << e.what()<<endl
              << "This message is normal and should be ignored"<<endl;
          }
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
