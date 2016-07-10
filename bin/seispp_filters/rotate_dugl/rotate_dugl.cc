/* This is a special purpose program to rotate 3C data for the dugl
   experiment at homestake.   Uses a straight line ray approximation
   that would normally be a dumb idea.  Also has hard wired station
   list of surface stations that are handled differently.  That is underground
   sites rotate to P aligned to straight ray, radial perpendicular in the
   ray vertical plan, and transverse the horizontal perpendicular to that
   plane.  For surface sites only the horizontals are rotate. */
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
    cerr << "rotate_dugl < infile > outfile"
        <<endl
        << "infile and outfile are a single ThreeComponentEnsemble boost serialization file"<<endl
        << "special purpose 3C rotation to P,R,T for DUGL (homestake) experiment"<<endl
        << "Use this only as a starting point for modification for any other data"
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
set<string> surface_sites()
{
    set<string> ss;
    ss.insert("YATES");
    ss.insert("WTP");
    ss.insert("ORO");
    ss.insert("ROSS");
    ss.insert("LHS");
    ss.insert("RRDG");
    return ss;
}
/* This obnoxious external variable is a necessary evil to deal with 
   error logging in the SEISPP library. Your code will probably not link 
   without it.*/ 
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    if(argc==1) usage();
    set<string> surfsta=surface_sites();
    set<string>::iterator ssptr;
    try{
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentEnsemble d;
        d=read_object<ThreeComponentEnsemble>(ia);
        vector<ThreeComponentSeismogram>::iterator dptr;
        int k;
        for(dptr=d.member.begin();dptr!=d.member.end();++dptr)
        {
            string sta=dptr->get_string("sta");
            double r[3],s[3],nu[3];
            double phi;
            r[0]=dptr->get_double("rx");
            r[1]=dptr->get_double("ry");
            r[2]=dptr->get_double("rz");
            s[0]=dptr->get_double("sx");
            s[1]=dptr->get_double("sy");
            s[2]=dptr->get_double("sz");
            double offset(0.0);
            for(k=0;k<3;++k)
            {
                double dx=r[k]-s[k];
                nu[k]=dx;
                offset += dx*dx;
            }
            offset=sqrt(offset);
            for(k=0;k<3;++k) nu[k]=nu[k]/offset;
            ssptr=surfsta.find(sta);
            if(ssptr==surfsta.end())
            {
                phi=atan2(nu[0],nu[1]);
                dptr->rotate(phi);
            }
            else
                dptr->rotate(nu);
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


