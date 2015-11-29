#include <iostream>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include "stock.h"
#include "perf.h"
#include "seispp.h"
#include "LatLong-UTMconversion.h"
/* This could be in an include, but will insert this prototype here
   rather than make an include file with one line.*/
TimeSeries ReadSegyTrace(FILE *);

void usage()
{
    cerr << "SU3CEnsembleConverter outfile [-pf pffile] < SUfile "<<endl
        << "Translates one ensemble of 3C data (sensor order 1,2,3=sensor 1)"
        <<endl
        << " in SU format to SEISPP ThreeComponentEnsemble object serialization"
        <<endl
        << " SU data comes through stdin - normally pipeline with an suwind filer"
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
    if(argc<2)usage();
    xargc=argc;
    xargv=argv;
    string outfile(argv[1]);
    string pffile("SU3CEnsembleConverter");
    int i,j;
    for(i=2;i<argc;++i)
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
        /* Now build all the list of metadata to be loaded */
        MetadataList tmdl=pfget_mdlist(pf,"trace_metadata_list");
        Metadata control(pf);
        bool apply_spreading_correction=control.get_bool("correct_for_geometric_spreading");
        double spow(0.0);  
        if(apply_spreading_correction)
        {
            spow=control.get_double("spreading_power_factor");
            cout << "Applying geometric spreading correction using "
                << "offset^"<<spow<<endl;
        }
        bool apply_rotation=control.get_bool("apply_rotation");
        double rotation_angle(0.0);
        if(apply_rotation)
        {
            rotation_angle=control.get_double("rotation_angle");
            /* All angles in my library are radians but input is degrees*/
            cout << "Horizontal components will be rotated by phi="<<rotation_angle<<endl;
            rotation_angle=rad(rotation_angle);
        }
        bool convert_utm=control.get_bool("convert_utm_to_dd");
        string zone;  //utm zone when needed
        int RefEl(23);   //frozen as WGS-84
        if(convert_utm)
        {
            zone=control.get_string("UTM_zone");
            cout << "UTM conversion will be done with zone="<<zone<<endl;
        }
        else
            cout << "Warning:  utm conversion is off."<<endl
                << "Using the result in ParticleMotionVTKconverter will fail"<<endl;
        /* This is a crude way to define the hang,vang orientation data. Works
        only for homestake data where all sensors had a common orientation */
        Tbl *t;
        t=pfget_tbl(pf,const_cast<char *>("channel_orientation"));
        if(t==NULL)
        {
            cerr << "Missing required Tbl parameter channel_orientation"<<endl;
            usage();
        }
        if(maxtbl(t)!=3)
        {
            cerr << "channel_orientation Tbl parameter in pf file is incorrectly defined"
                <<endl << "Size must be 3.   Size found="<<maxtbl(t)<<endl;
            exit(-1);
        }
        double hang[3],vang[3];   // parallel arrays
        long int ii;
        for(ii=0;ii<3;++ii) 
        {
            char *line;
            line=(char *)gettbl(t,ii);
            sscanf(line,"%d%lf%lf",&i,hang+ii,vang+ii);
            if(i!=ii)
            {
                cerr << "Format error for channel_orientation data"<<endl
                    << "Must be listed in ascending channel order (0,1,2)"<<endl
                    << "Read index "<<i<<" when expecting "<<ii<<endl;
                exit(-1);
            }
            cout << "Component "<<ii<<" Setting hang="<<hang[ii]<< " and vang="
                <<vang[ii]<<endl;
        }

        /* This will hold our results */
        ThreeComponentEnsemble ens;
        cout << "SU3CEnsembleConverter processing begins - reading from stdin"<<endl;
        /* Now load the data file until EOF.  This is not general
           but assumes this program will always follow a suwind 
           command to build the ensemble.  The logic of this
           read loop is a bit perverted because of the read procedure
           interface to SU was designed to return a TimeSeries object.
           */
        bool readok;
        TimeSeries dread=ReadSegyTrace(stdin);
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
            dread.put("hang",hang[k]);
            dread.put("vang",vang[k]);
            channels[k]=dread;
            if(k==2)
            {
                ThreeComponentSeismogram d3c(channels,0);
                if(apply_rotation)
                    d3c.rotate(rotation_angle);
                if(convert_utm)
                {
                    double easting,northing;
                    easting=d3c.get_double("rx");
                    northing=d3c.get_double("ry");
                    double lat,lon;
                    /* This library routine returns lat,lon in degrees.
                       We post that to header as css names to mate
                       with ParticleMotionVTKconverter.  */
                    UTMtoLL(RefEl,northing,easting,zone.c_str(),lat,lon);
                    d3c.put("site.lat",lat);
                    d3c.put("site.lon",lon);
                }
                if(apply_spreading_correction)
                {
                    double offset=d3c.get_double("offset");
                    offset=fabs(offset);
                    /* Test for small offset - unit dependent assumes
                       offset in m */
                    if(offset>0.01)
                    {
                        double scale=pow(offset,spow);
                        dscal(3*d3c.ns,scale,d3c.u.get_address(0,0),1);
                    }
                }
                ens.member.push_back(d3c);
            }
            dread=ReadSegyTrace(stdin);
            if(dread.ns<=0) readok=false;
            ++n;
            k=n%3;
        }
        cout << "Finished bundling 3C ensemble with "<<ens.member.size()
            << " 3C objects"<<endl
        << "Writing to archive file "<<outfile<<endl;
        /* Now we save the result to outfile */
        oa << ens;
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
}
