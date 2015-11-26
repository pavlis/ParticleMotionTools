#include <iostream>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include "stock.h"
#include "seispp.h"
#include "HeaderMap.h"
#include "AttributeCrossReference.h"
/* This could be in an include, but will insert this prototype here.
   In a separate file as this routine may eventually become a library
   routine */
TimeSeries ReadSegyTrace(FILE *, HeaderMap& hm, AttributeCrossReference& xref,
        MetadataList& mdl);

/* This procedure parses a Tbl with the tag "SUAttributeCrossReference" 
   and constructs the cross reference map needed by the GenericFileHandle. 
   Useful here to allow names that are more rational for passing to 
   output object.   

   Input is converted to a large string and passed to the constructor.
   Returns the constructed object.   Can throw an exception so 
   there is a catch all handler. 
   */
AttributeCrossReference parse_xref_tbl(Pf *pf)
{
    try{
        Tbl *t=pfget_tbl(pf,const_cast<char *>("SUAttributeCrossReference"));
        if(t==NULL)
            throw SeisppError(string("parse_xref_tbl:   pf is missing required")
                    +" Tbl with tab SUAttributeCrossReference");
        string xrefbuf;
        char *s;
        int i;
        for(i=0;i<maxtbl(t);++i)
        {
            s=(char *)gettbl(t,i);
            if(i==0)
                xrefbuf=string(s);
            else
                xrefbuf+=string(s);
            xrefbuf+="\n";
        }
        return AttributeCrossReference(xrefbuf);
    }catch(...){throw;};
}
void usage()
{
    cerr << "SU3CEnsembleConverter SUinfile outfile [-pf pffile]"<<endl
        << "Translates one ensemble of 3C data (sensor order 1,2,3=sensor 1)"
        <<endl
        << " in SU format to SEISPP ThreeComponentEnsemble object serialization"
        <<endl;
    exit(-1);
}
/* This is defined in SU's par.h and seems necessary for this to link.
      Not used in this code, but an annoying extern.   This is 
      an incredibly obscure trick to make this work.  Found by 
     pure hacking
*/
extern "C"{
int xargc;
char **xargv;
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    if(argc<3)usage();
    xargc=argc;
    xargv=argv;
    string SUinfile(argv[1]);
    string outfile(argv[2]);
    string pffile("SU3CEnsembleConverter");
    int i,j;
    for(i=3;i<argc;++i)
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
    Pf *pf;
    if(pfread(const_cast<char*>(pffile.c_str()),&pf))
    {
        cerr << "pfread failed on pffile="<<pffile<<endl;
        usage();
    }
    try{
        /* Make sure the output file can be openned immediately */
        std::ofstream ofp(outfile.c_str(),ios::out);
        if(ofp.fail()) 
        {
            cerr << "FATAL:  Open failed for file="<<outfile<<endl;
            exit(-1);
        }
        boost::archive::text_oarchive oa(ofp);
        /* This beast provides a cross reference between su names
           and SEISPP name conventions needed downstream */
        AttributeCrossReference xref=parse_xref_tbl(pf);
        /* Now build all the list of metadata to be loaded */
        MetadataList tmdl=pfget_mdlist(pf,"trace_metadata_list");
        Metadata control(pf);
        bool apply_rotation=control.get_bool("apply_rotation");
        double rotation_angle(0.0);
        if(apply_rotation)
        {
            rotation_angle=control.get_double("rotation_angle");
            /* All angles in my library are radians but input is degrees*/
            rotation_angle=rad(rotation_angle);
        }
        HeaderMap hm(pf,string("SEGYfloat"));
        /* This will hold our results */
        ThreeComponentEnsemble ens;
        /* Now load the data file until EOF.  This is not general
           but assumes this program will always follow a suwind 
           command to build the ensemble.  The logic of this
           read loop is a bit perverted because of the read procedure
           interface to SU was designed to return a TimeSeries object.
           */
        bool readok;
        TimeSeries dread=ReadSegyTrace(stdin,hm,xref,tmdl);
        if(dread.ns<=0) 
        {
            cerr << "No data to process.   ReadSegyTrace hit EOF immediately"
                <<endl;
            exit(-1);
        }
        else
            readok=true;
        int n;   // counter for number of traces read - sanity check
        int k;   // mod 3 to get channel code
        /* This holds channels to build 3c objects.  STL container
        has to be initialized like this to allow use of indexing
        operator */
        vector<TimeSeries> channels;
        for(n=0;n<3;++n) channels.push_back(TimeSeries());
        n=0;  k=0;
        while(readok)
        {
            double dcoord;
            int icoord;
            channels[k]=dread;
            if(k==2)
            {
                ThreeComponentSeismogram d3c(channels,0);
                if(apply_rotation)
                    d3c.rotate(rotation_angle);
                /* This is necessary because of the obnoxious
                   coordinate scaling in segy an storage of
                   coordinates as int.  Here we translate them
                   to floats and post the result to a new name */
                int scalco=d3c.get_int("scalco");
                double scale=(double)scalco;
                icoord=d3c.get_int("rx");
                dcoord=((double)icoord)/scale;
                d3c.put("rx_flaot",dcoord);
                icoord=d3c.get_int("ry");
                dcoord=((double)icoord)/scale;
                d3c.put("ry_flaot",dcoord);
                icoord=d3c.get_int("sx");
                dcoord=((double)icoord)/scale;
                d3c.put("sx_flaot",dcoord);
                icoord=d3c.get_int("sy");
                dcoord=((double)icoord)/scale;
                d3c.put("sy_flaot",dcoord);
                ens.member.push_back(d3c);
                k=0;
            }
            dread=ReadSegyTrace(stdin,hm,xref,tmdl);
            if(dread.ns<=0) readok=false;
        }
        /* Now we save the result to outfile */
        oa << ens;
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
}
