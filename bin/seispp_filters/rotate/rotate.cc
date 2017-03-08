/* General purpose routine to rotate an ensemble.   Works in one
   of two modes depending on switch constant_transformation.  When true
   rotates to a specified orientation defined by spherical angles phi and
   theta (see ThreeComponentSeismogram rotate method for details).  When
   false the program will attempt to extract source and receiver coordinates
   from each seismogram to compute backazimuth.  Then the horizontals are
   rotated to radial-transverse.   

   Does not currently support something like a P wave emergence angle formula
   assuming that is better done with the free-surface transformation method.
   May want to actually implement that as an option. */
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
void usage()
{
    cerr << "rotate_dugl [-pf pffile] < infile > outfile"
        <<endl
        << "infile is one the concatenation of one or more ThreeComponentEnsemble objects"
        <<endl
        << "Use -pf to specify alternate parameter file to default rotate.pf"
        <<endl;
    exit(-1);
}
class RotateControl
{
  public:
    bool constant_transformation;
    double phi,theta;
    string stalatkey,stalonkey;
    string evlatkey,evlonkey;
    RotateControl(string pffile);
};
RotateControl::RotateControl(string pffile)
{
  Pf *pf;
  if(pfread(const_cast<char *>(pffile.c_str()),&pf))
    throw SeisppError("RotateControl:  pfread failed for "+pffile);
  try{
    Metadata md(pf);
    constant_transformation=md.get_bool("constant_transformation");
    if(constant_transformation)
    {
      phi=md.get_double("phi");
      theta=md.get_double("theta");
      /* convert both to radians */
      phi=rad(phi);
      theta=rad(theta);
    }
    else
    {
      phi=0.0;
      theta=0.0;
      stalatkey=md.get_string("station_latitude_key");
      stalonkey=md.get_string("station_longitude_key");
      evlatkey=md.get_string("event_latitude_key");
      evlonkey=md.get_string("event_longitude_key");
    }
  }catch(...){throw;};
}

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
  string pffile("rotate");
  for(i=1;i<argc;++i)
  {
    string sarg(argv[i]);
    if(sarg=="-pf")
    {
      ++i;
      if(i>=argc) usage();
      pffile=string(argv[i]);
    }
    else
      usage();
  }
    try{
        RotateControl rc(pffile);
        SphericalCoordinate ctsc;
        ctsc.radius=1.0;
        ctsc.phi=rc.phi;
        ctsc.theta=rc.theta;
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentEnsemble d;
        d=read_object<ThreeComponentEnsemble>(ia);
        vector<ThreeComponentSeismogram>::iterator dptr;
        int k;
        for(dptr=d.member.begin();dptr!=d.member.end();++dptr)
        {
          double slat,slon,evlat,evlon;
          double az,delta;
          if(rc.constant_transformation)
          {
            dptr->rotate(ctsc);
          }
          else
          {
            slat=dptr->get_double(rc.stalatkey);
            slon=dptr->get_double(rc.stalonkey);
            evlat=dptr->get_double(rc.evlatkey);
            evlon=dptr->get_double(rc.evlonkey);
            dist(slat,slon,evlat,evlon,&delta,&az);
            dptr->rotate(az+M_PI_2);
          }
        }
        write_object<ThreeComponentEnsemble>(d,oa);
    }catch(boost::archive::archive_exception const& e)
    {
        cerr << "Error in read or write of serialization object"<<endl
            << "This was the message posted:"<<endl
            << e.what()<<endl;
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(std::exception& stexc)
    {
       cerr << stexc.what()<<endl;
    }
}
