#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "stock.h"
#include "seispp.h"
#include "HFArray.h"
#include "PMTimeSeries.h"
#include "TimeWindow.h"
#include "PfStyleMetadata.h"
#include "AttributeMap.h"
using namespace std;
using namespace SEISPP;
vector<PMTimeSeries> read_pmdata(string listfile)
{
        vector<PMTimeSeries> result;
        ifstream lfin;
        lfin.open(listfile.c_str());
        if(!lfin.good())
        {
            cerr << "Open failed on input list file="<<listfile<<endl
                << "Cannot recover - exiting"<<endl;
            exit(-1);
        }
        char fname[1024];
        while(lfin.getline(fname,1024))
        {
            ifstream ifs(fname);
            if(ifs.good())
            {
                boost::archive::text_iarchive ar(ifs);
                PMTimeSeries pmtsin;
                try {
                    ar >> pmtsin;
                    result.push_back(pmtsin);
                }catch(std::exception& err)
                {
                    cerr << "Error reading archive file = "<<fname<<endl
                        << "Message posted by boost::archive:  "
                        <<err.what()<<endl;
                }
                lfin.close();
                ifs.close();
            }
            else
            {
                cerr << "Open failure for input archive file="
                    << fname <<endl
                    << "File skipped.   Trying net on list "<<endl;
            }
        }
        lfin.close();
        return result;
}
        
/* This procedure is in vtp_writer.cc*/
void WriteEllipses(vector<PMTimeSeries>& d,
        double t, double scale, int np_per_ellipse,
                ofstream& out);

void write_animation_files(vector<PMTimeSeries>& pmdv,
        TimeWindow& dtw,double scale, int npoints_ellipse,
        string dir, string fbase, string pvdfile)
{
    const string base_error("write_animation_files:  ");
    const string fext(".vtp");
    try {
        if(makedir(const_cast<char *>(dir.c_str())))
            throw SeisppError(base_error 
                    + "makedir failed to create directory="
                    + dir);
        int filenumber(1);  // file name tag - incremented for each 
        list<string> allfnames;
        string path;
        ofstream fout;
        /* assume first member is representative.   Caller
           should assure input is rational and start times and
           data lengths are equal or at least close.  */
        int number_time_steps=trunc(dtw.length()/pmdv[0].dt);
        double dt=pmdv[0].dt;
        cout << "Program is writing data for "<< number_time_steps 
            << " time steps at time intervals of "
            << dt<<endl;
        int i;
        double t;
        for(int i=0,t=dtw.start;i<number_time_steps;++i)
        {
            stringstream ss;
            ss << dir <<"/"<<fbase<<"_"<<filenumber<<fext;
            string fname=ss.str();
            fout.open(fname.c_str(),std::ofstream::out);
            if(fout.fail())
            {
                cerr << base_error << "open failed for file="
                    <<fname<<endl
                    << "No output for time step "<<i<<endl;
            }
            else
            {
                try{
                    WriteEllipses(pmdv,dtw.start,scale,npoints_ellipse,fout);
                }catch(SeisppError& serr)
                {
                    cerr << "WriteTimeWindow threw this error for time step "
                        << i<<endl;
                    serr.log_error();
                    continue;
                }
                allfnames.push_back(fname);
                ++filenumber;
                dtw=dtw.shift(dt);
                fout.close();
            }
        }
        /* Now write the pvd file */
        fout.open(pvdfile.c_str(),std::ofstream::out);
        if(fout.fail())
        {
            cerr << "open failed for pvd file="<<pvdfile<<endl;
            cerr << "You need to build this by hand or rerun"<<endl;
            exit(-1);
        }
        /* preamble */
        fout << "<?xml version=\"1.0\"?>"<<endl 
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"
            <<endl
            << "<Collection>"<<endl;
        
        list<string>::iterator iptr;
        for(i=0,iptr=allfnames.begin();iptr!=allfnames.end();++iptr,++i)
        {
            fout << "<DataSet timestep=\"" << i <<"\" file=\""
                << *iptr <<"\"/>" << endl;
        }
        fout << "</Collection>"<<endl
            << "</VTKFile>"<<endl;
        fout.close();
    }catch(...){throw;};
}
double find_max_amplitude(PMTimeSeries& d)
{
    int i;
    double maxamp(0.0);
    for(i=0;i<d.ns;++i)
    {
        ParticleMotionEllipse e=d.ellipse(i);
        if(e.majornrm>maxamp) maxamp=e.majornrm;
    }
    return maxamp;
}


void show_scaling(vector<PMTimeSeries>& pmdv,
        double initial_scale_factor)
{
    vector<double> ranges,low,high,work;
    double ampmax,awork;
    int i;
    vector<PMTimeSeries>::iterator dptr;
    for(dptr=pmdv.begin(),i=0;dptr!=pmdv.end();++dptr,++i)
    {
        if(i==0)
        {
            low.push_back(dptr->get_double("x"));
            low.push_back(dptr->get_double("y"));
            low.push_back(dptr->get_double("z"));
            high=low;
            work=low;
            ampmax=find_max_amplitude(*dptr);
        }
        else
        {
            work[0]=dptr->get_double("x");
            work[1]=dptr->get_double("y");
            work[2]=dptr->get_double("z");
            for(i=0;i<3;++i)
            {
                if(work[i]>high[i]) high[i]=work[i];
                if(work[i]<low[i]) low[i]=work[i];
            }
            awork=find_max_amplitude(*dptr);
            if(awork>ampmax) ampmax=awork;
        }
    }
    for(i=0;i<3;++i) 
        ranges.push_back(high[i]-low[i]);
    cout << "PMTimeSeriesToVTK:   Current scale factor="
        << initial_scale_factor <<endl
        << " produces a maximum, scaled particle motion amplitude = "
        << ampmax*initial_scale_factor << endl
        << "Here are coordinate ranges (low, high, extent):"  << endl;
    for(i=0;i<3;++i)
        cout << "x"<<i+1<<":  "<<low[i]<<", "<<high[i]<<", "<<ranges[i]<<endl;
}
void load_cartesian(vector<PMTimeSeries>& d, RegionalCoordinates& coords)
{
    try {
        vector<PMTimeSeries>::iterator dptr;
        for(dptr=d.begin();dptr!=d.end();++dptr)
        {
            double lat,lon,elev,r;
            lat=dptr->get_double("site.lat");
            lon=dptr->get_double("site.lon");
            elev=dptr->get_double("site.elev");
            /* Because of site qualifier assume these degree units */
            lat=rad(lat);
            lon=rad(lon);
            r=r0_ellipse(lat)+elev;
            Cartesian_point cp=coords.cartesian(lat,lon,r);
            dptr->put("x",cp.x1);
            dptr->put("y",cp.x2);
            dptr->put("z",cp.x3);
        }
    }catch(...){throw;};
}


void usage()
{
    cerr << "PMTimeSeriesToVTK  inlist [-pf pffile]"<<endl
        << "inlist is list of serialized file names for PMTimeSeries object"
        <<endl
        << "Use -pf to alter default control file PMTimeSeriesToVTK.pf"
        <<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc<2) usage();
    string infilelist(argv[1]);
    string pffile("PMTimeSeriesToVTK.pf");
    int i;
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
    try {
        PfStyleMetadata md=pfread(pffile);
        bool sap_mode=md.get_bool("small_aperture_array_mode");
        HFArray *array;
        if(sap_mode)
        {
            string array_file=md.get_string("array_geometry_file");
            array=new HFArray(array_file);
        }
        MetadataList tracemdl=get_mdlist(md,"trace_mdlist");
        MetadataList ensemblemdl=get_mdlist(md,"ensemble_mdlist");
        AttributeMap am("css3.0");
        /* This defines the relative time window , which must 
           be shorter or equal to data available */
        TimeWindow dtw;
        dtw.start=md.get_double("display_window_start_time");
        dtw.end=md.get_double("display_window_end_time");
        /* This defines file names for output files.   
           animate_file_base - base file name for animation pieces
                (These get number tags for each frame)
           animate_file_directory - output directory for animation files
           pvd_file - wrapper file that defines how animate_file_base files
                are to be combined for animation */
        string animate_file_base=md.get_string("animate_file_base");
        string pvd_file=md.get_string("pvd_file");
        string animate_file_directory=md.get_string("animate_file_directory");
        int npoints=md.get_int("number_points_per_ellipse");
        /* This builds the object used to convert data to a local
           coordinate system for paraview */
        double olat,olon,odepth,azn;
        olat=md.get_double("origin_latitude");
        olon=md.get_double("origin_longitude");
        odepth=md.get_double("origin_depth");
        azn=md.get_double("origin_azimuth_north");
        olat=rad(olat);  olon=rad(olon); azn=rad(azn);
        RegionalCoordinates coords(olat,olon,r0_ellipse(olat)-odepth,azn);
        /* This is the scale factor used to convert raw sample data
           to be visible figures at map scale. Note interactive 
           loop below used because this is a critical parameter */
        double PMscale_factor=md.get_double("PMscale_factor");
        /* Read the data defined by the input list of files */
        vector<PMTimeSeries>::iterator dptr;
        vector<PMTimeSeries> d=read_pmdata(infilelist);
        if(sap_mode)
        {
            /* In small aperture array mode we override any
               receiver position information with coordinate data
               stored in the HFArray object */
            for(dptr=d.begin();
                    dptr!=d.end();++dptr)
            {
                try {
                    string sta=dptr->get_string("sta");
                    Geographic_point gp=array->geographic_location(sta);
                    /* Metadata coordinates are in degrees, but 
                       gp has units of radians */
                    dptr->put("site.lat",deg(gp.lat));
                    dptr->put("site.lon",deg(gp.lon));
                    /* elev has to be computed from radius like this.
                       This is repeated elsewhere, which is
                       inefficient, but necessary to have this work for
                       both small aperture and large arrays */
                    double elev;
                    elev=gp.r-r0_ellipse(gp.lat);
                    dptr->put("site.elev",elev);
                }catch(SeisppError& serr)
                {
                    cerr << "Failure in loop setting coordinates in "
                        <<"small aperture array mode."
                        <<"SeisppError message posted:"<<endl;
                    serr.log_error();
                    cerr << "Additional errors are likely"<<endl;
                }
            }
        }
        /* This procedure loads metadata x,y,z attributes
           as Cartesian (km) coordinates using site.lat, site.lon,
           and site.elev.   */
        load_cartesian(d,coords);
        double newscale;
        newscale=PMscale_factor;
        /* Scan the data for largest amplitude and set the scale
           factor interactively. */
        string ques;
        do {
            show_scaling(d,newscale);
            cout << "Answer y to accept this scale "
                    << "(anything else to try again): "
                    <<endl;
            cin >> ques;
            if(ques!="y")
            {
                cout << "Enter new scale factor->";
                cin >> newscale;
            }
        }while(ques!="y");
        /* Now we write the string of files for animation.   
           This procedure does that */
        write_animation_files(d,dtw,newscale,npoints,
                animate_file_directory,animate_file_base,pvd_file);
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(...)
    {
        cerr << "Something threw an unexpected exception"
            << "Fatal - exiting"
            <<endl;
    }
}
