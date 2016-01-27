/* This set of system includes are always required.  Do not remove them.*/
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
/* These are SEISPP includes that are a base requirement 
   to use this template file */
#include "seispp.h"
/* Replace these as needed.  This one is here only to make
   this do nothing template compile for testing configuration */
#include "ensemble.h"
#include "PfStyleMetadata.h"
/* You will get lots of errors without these namespace
   declaration*/
using namespace std;   // most compilers do not require this
using namespace SEISPP;  //This is essential to use SEISPP library
class EnsembleCoordinates
{
public:
    EnsembleCoordinates(PfStyleMetadata& control);
    void set_coords(ThreeComponentEnsemble& d);
private:
    vector<double> rx,ry,sx,sy,relev,selev,offset;
};
EnsembleCoordinates::EnsembleCoordinates(PfStyleMetadata& md)
{
  try {
    list<string> lines;
    lines=md.get_tbl("coordinates");
    list<string>::iterator lptr;
    int i;
    for(i=0,lptr=lines.begin();lptr!=lines.end();++i,++lptr)
    {
        /* Assuming each line has order:  rx ry relev sx sy selev*/
        stringstream ss(*lptr);
        double d;
        ss >>d;
        rx.push_back(d);
        ss >>d;
        ry.push_back(d);
        ss >>d;
        relev.push_back(d);
        ss >>d;
        sx.push_back(d);
        ss >>d;
        sy.push_back(d);
        ss >>d;
        selev.push_back(d);
        d=(rx[i]-sx[i])*(rx[i]-sx[i]) + (ry[i]-sy[i])*(ry[i]-sy[i]);
        offset.push_back(sqrt(d));
    }
  }catch(...){throw;};
}
void EnsembleCoordinates::set_coords(ThreeComponentEnsemble& d)
{
    /* All the vectors have to be the same length so we only need
       test one for a sanity check */
    if(offset.size()!=d.member.size())
    {
        stringstream sserr;
        sserr << "EnsembleCoordinates::set_coords:   size mismatch"
            << endl
            << "Coordinate table size="<<offset.size()
            << " but ensemble size = "<< d.member.size()<<endl;
        throw SeisppError(sserr.str());
    }
    int i;
    vector<ThreeComponentSeismogram>::iterator dptr;
    if(SEISPP_verbose) cerr << "member rx   ry   relev  sx   sy   selev   offset"<<endl;
    for(i=0,dptr=d.member.begin();dptr!=d.member.end();++dptr,++i)
    {
        dptr->put("rx",rx[i]);
        dptr->put("ry",ry[i]);
        dptr->put("relev",relev[i]);
        dptr->put("sx",sx[i]);
        dptr->put("sy",sy[i]);
        dptr->put("selev",selev[i]);
        dptr->put("offset",offset[i]);
        if(SEISPP_verbose)
        {
            cerr << i <<" "
                << dptr->get_double("rx") << " "
                << dptr->get_double("ry") << " "
                << dptr->get_double("relev") << " "
                << dptr->get_double("sx") << " "
                << dptr->get_double("sy") << " "
                << dptr->get_double("selev") << " "
                << dptr->get_double("offset") << endl;
        }
    }
}

void usage()
{
    cerr << "set_coords [-pf pffile] < infile > outfile"
        <<endl
        << "sets coordinates in a 3C ensemble driven by pffile table"<<endl
        << "Table must be in channel order grouped by threes - very crude"
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
    int i;
    const int narg_required(0);
    string pffile("set_coords.pf");
    for(i=narg_required+1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-pf")
        {
            ++i;
            if(i>=argc)usage();
            pffile=string(argv[i]);
        }
        else
            usage();
    }
    try{
        PfStyleMetadata control=pfread(pffile);
        EnsembleCoordinates coords(control);
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
        coords.set_coords(d);
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

