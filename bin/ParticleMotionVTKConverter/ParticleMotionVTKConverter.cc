#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "stock.h"
#include "seispp.h"
#include "ParticleMotionData.h"
#include "dbpp.h"
#include "ensemble.h"
#include "StationChannelMap.h"
#include "filter++.h"
#include "HFArray.h"
/* This procedure is in another file and this should probably be an include */
void WriteTimeWindow(vector<ParticleMotionData>& d,
        TimeWindow tw,ofstream& out);
using namespace std;
using namespace SEISPP;
ThreeComponentEnsemble *get_event_data(DatascopeHandle& dbh,
        long evid, MetadataList& emdl, MetadataList& tmdl,
        AttributeMap& am)
{
    try{
        /* Dogmatically use wfprocess for this read approach */
        dbh.lookup("evlink");
        stringstream ss;
        ss << "evid == "<<evid;
        dbh.subset(ss.str());
        dbh.natural_join("wfprocess");
        dbh.natural_join("event");
        dbh.natural_join("origin");
        dbh.subset(string("orid==prefor"));
        dbh.natural_join("assoc");
        dbh.natural_join("arrival");
        dbh.natural_join("site");
        int nmembers=dbh.number_tuples();
        ThreeComponentEnsemble *d;
        d=new ThreeComponentEnsemble();
        d->member.reserve(nmembers);
        /* this is ugly - copied from database constructor.
        Point is to create ensemble metadata */
        Metadata ensmd(dbh,emdl,am);
        copy_selected_metadata(ensmd,dynamic_cast<Metadata&>(*d),emdl);
        /* Now we actually read the data */
        for(int i=0;i<nmembers;++i,++dbh)
        {
            ThreeComponentSeismogram *s;
            try {
                s=new ThreeComponentSeismogram(dbh,tmdl,am);
            }catch(SeisppError& serr)
            {
                serr.log_error();
                cerr << "Skipping data for that station"<<endl;
                delete s;
                continue;
            }
            d->member.push_back(*s);
            delete s;
        }
        return d;
    }catch(...){throw;};
}
void set_ensemble_scale(vector<ParticleMotionData>& pmdv,double scale)
{
    vector<ParticleMotionData>::iterator dptr;
    for(dptr=pmdv.begin();dptr!=pmdv.end();++dptr)
        dptr->set_scale(scale);
}
/* This procedure computes and writes to cout a set of parameters
   that can be used for an interactive decision on scaling.   
   Note there are no error handlers here because the current
   implementation does not throw errors when requesting data on
   the full data extent.   If altered to a different time window 
   this could be an issue.
*/ 
void show_scaling_parameters(vector<ParticleMotionData>& pmdv)
{
    vector<double> ranges,low,high,work;
    double ampmax,awork;
    int i;
    vector<ParticleMotionData>::iterator dptr;
    for(dptr=pmdv.begin(),i=0;dptr!=pmdv.end();++dptr,++i)
    {
        if(i==0)
        {
            low=dptr->center();
            high=low;
            ampmax=dptr->max_amplitude();
        }
        else
        {
            work=dptr->center();
            /* Assume work size is 3 */
            for(i=0;i<3;++i)
            {
                if(work[i]>high[i]) high[i]=work[i];
                if(work[i]<low[i]) low[i]=work[i];
            }
            awork=dptr->max_amplitude();
            if(awork>ampmax) ampmax=awork;
        }
    }
    for(i=0;i<3;++i) 
        ranges.push_back(high[i]-low[i]);
    cout << "ParticleMotionVTKConverter:   Current scale factor="
        << pmdv[0].get_scale()<<endl
        << "This produces a maximum, scaled particle motion amplitude = "
        << ampmax << endl
        << "Here are coordinate ranges (low, high, extent):"  << endl;
    for(i=0;i<3;++i)
        cout << "x"<<i+1<<":  "<<low[i]<<", "<<high[i]<<", "<<ranges[i]<<endl;
}
void write_animation_files(vector<ParticleMotionData>& pmdv,
        double twlen, double dt, string dir, string fbase, string pvdfile,
        bool engine_mode)
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
        /* Assume all data have the same time span - may be dangerous */
        TimeWindow fulltimewindow=pmdv[0].timespan();
        double trange=fulltimewindow.length();
        /* This computes the number of time steps and asks for confirmation
           as a sanity check */
        int number_time_steps;
        string test;
        do {
            number_time_steps=trunc( (fulltimewindow.end-twlen)/dt );
            cout << "write_animation_files procedure:   "
                << "Data time range is "
                   <<fulltimewindow.start<<" to "<<fulltimewindow.end<<endl
                << "Current animation time window length="<<twlen
                << " with animation time step="<<dt<<endl
                << "With these the computed number of time steps is "
                << number_time_steps<<endl;
            if(engine_mode)
            {
                test="y";
            }
            else
            {
                cout << "Ok to proceed writing files "
                    <<"(Enter y to accept - anything means no):";
                cin >> test;
                if(test!="y")
                {
                    cout << "Enter new value for animation time window length: ";
                    cin >> twlen;
                    cout << "Enter new value for animation time step: ";
                    cin >> dt;
                }
            }
        }while(test!="y");
        /* This is the moving time window used for extraction */
        TimeWindow tw(fulltimewindow.start,fulltimewindow.start+twlen);
        for(int i=0;i<number_time_steps;++i)
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
                WriteTimeWindow(pmdv,tw,fout);
                }catch(SeisppError& serr)
                {
                    cerr << "WriteTimeWindow threw this error for time step "
                        << i<<endl;
                    serr.log_error();
                    continue;
                }
                allfnames.push_back(fname);
                ++filenumber;
                tw=tw.shift(dt);
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
        int i;
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


void usage()
{
    cerr << "ParticleMotionVTKConverter db "
        <<" (-a start_time end_time | -e evid) [-filter fstring -engine"
        <<" -pf pffile]"<<endl
        << "Must use either -a for fixed time window or -e for one event"<<endl
        << "-filter overrides pf definition of default filter"<<endl
        << "-engine disables interactive scale editing"<<endl
        << "Use -pf to alter default control file ParticleMotionVTKConverter.pf"
        <<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc<2) usage();
    Pf *pf;
    if(pfread(const_cast<char *>("ParticleMotionVTKConverter"),&pf))
    {
        cerr << "pfread failed trying to read default pf file="
            << "ParticleMotionVTKConverter.pf"<<endl
            << "This file should have been installed in $ANTELOPE/data/pf"
            <<endl;;
        usage();
    }
    string dbname(argv[1]);
    bool abstime(false),eventmode(false);
    bool engine_mode(false);
    TimeWindow twtotal;
    string filtername_from_arglist("none");
    long evid;
    int i;
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-a")
        {
            abstime=true;
            if(i>argc-2) 
                usage();
            ++i;
            twtotal.start = str2epoch(argv[i]);
            ++i;
            twtotal.end = str2epoch(argv[i]);
        }
        else if(sarg=="-e")
        {
            eventmode=true;
            ++i;
            if(i>=argc) usage();
            evid=atol(argv[i]);
        }
        else if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pffree(pf);
            if(pfread(argv[i]),&pf)
            {
                cerr << "pfread error using pffile="<<argv[i]<<endl;
                usage();
            }
        }
        else if(sarg=="-engine")
            engine_mode=true;
        else if(sarg=="-filter")
        {
            ++i;
            if(i>=argc) usage();
            sarg=string(argv[i]);
            filtername_from_arglist=sarg;
        }
        else
            usage();
    }
    if(!(eventmode || abstime) ) usage();
    try {
        Metadata md(pf);
        bool sap_mode=md.get_bool("small_aperture_array_mode");
        HFArray *array;
        if(sap_mode)
        {
            string array_file=md.get_string("array_geometry_file");
            array=new HFArray(array_file);
        }
        MetadataList tracemdl=pfget_mdlist(pf,"trace_mdlist");
        MetadataList ensemblemdl=pfget_mdlist(pf,"ensemble_mdlist");
        AttributeMap am("css3.0");
        /* This defines the relative time window saved, which must 
           be shorter than data available */
        TimeWindow tw_relative;
        tw_relative.start=md.get_double("process_window_start_time");
        tw_relative.end=md.get_double("process_window_end_time");
        /* This defines file names for output files.   
           full_window_file - used to write all data in process_window
           animate_file_base - base file name for animation pieces
                (These get number tags for each frame)
           animate_file_directory - output directory for animation files
           pvd_file - wrapper file that defines how animate_file_base files
                are to be combined for animation */
        string full_window_file=md.get_string("full_window_file");
        string animate_file_base=md.get_string("animate_file_base");
        string pvd_file=md.get_string("pvd_file");
        string animate_file_directory=md.get_string("animate_file_directory");
        /* These control animation output.  animate_twlength is the time
           window (in seconds) for each frame and animate_dt is the
           time step per frame. */
        double animate_dt=md.get_double("animate_dt");
        double animate_twlength=md.get_double("animate_time_window_length");
        bool filter_data=md.get_bool("filter_data");
        string filter_definition;
        if(filter_data)
            filter_definition=md.get_string("BRTT_filter_definition");
        /* override that name if arg list had a filter definition*/
        if(filtername_from_arglist!="none")
            filter_definition=filtername_from_arglist;
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

        ThreeComponentEnsemble *draw;
        auto_ptr<ThreeComponentEnsemble> d;
        double stime,etime;  /* Total time period to read */
        DatascopeHandle dbh(dbname,true);
        /* The prep routines form the working view with datascope.*/
        if(eventmode)
        {
            draw=get_event_data(dbh,evid,ensemblemdl,tracemdl,am);
        }
        else 
        {
            StationChannelMap scm(pf);
            draw=new ThreeComponentEnsemble(dbh,twtotal,scm);
        }
        /* Optional filter */
        vector<ThreeComponentSeismogram>::iterator dptr;
        if(filter_data)
        {
            TimeInvariantFilter f(filter_definition);
            for(dptr=draw->member.begin();
                    dptr!=draw->member.end();++dptr)
                f.apply(*dptr);
        }
        if(eventmode)
        {
            d=ArrivalTimeReference(*draw,string("arrival.time"),tw_relative);
        }
        else
        {
            d=auto_ptr<ThreeComponentEnsemble>
                (new ThreeComponentEnsemble(dynamic_cast<Metadata &>(*draw),
                                            draw->member.size()));
            for(dptr=draw->member.begin();
                    dptr!=draw->member.end();++dptr)
            {
                //DEBUG
                /*
                dmatrix utmp;
                utmp=tr(dptr->u);
                cout << "Raw Data "<<endl<< utmp <<endl;
                */
                ThreeComponentSeismogram work(*dptr);
                work.ator(dptr->t0);
                //DEBUG
                /*
                utmp=tr(work.u);
                cout << "After call to ator"<<endl<<utmp<<endl;
                */
                d->member.push_back(WindowData(work,tw_relative));
            }
        }
        delete draw;
        if(sap_mode)
        {
            /* In small aperture array mode we override any
               receiver position information with coordinate data
               stored in the HFArray object */
            for(dptr=d->member.begin();
                    dptr!=d->member.end();++dptr)
            {
                try {
                    string sta=dptr->get_string("sta");
                    Geographic_point gp=array->geographic_location(sta);
                    /* Metadata coordinates are in degrees, but 
                       gp has units of radians */
                    dptr->put("site.lat",deg(gp.lat));
                    dptr->put("site.lon",deg(gp.lon));
                    /* elev has to be computed from radius like this.
                       This is repeated in ParticleMotionData, which is
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
        /* Now we need to convert the ensemble to a vector of
           the ParticleMotionData objects that are used for 
           building vtp files*/
        vector<ParticleMotionData> pmdv;
        pmdv.reserve(d->member.size());
        for(dptr=d->member.begin();dptr!=d->member.end();++dptr)
            pmdv.push_back(ParticleMotionData(*dptr,coords));
        d.reset();
        double newscale;
        newscale=PMscale_factor;
        /* First write the file containing the full data window */
        string ques;
        do {
            ofstream outall;
            outall.open(full_window_file.c_str(),std::ofstream::out);
            if(outall.fail())
            {
                cerr << "Open failed on full_window_file="
                    << full_window_file<<endl
                    << "Fatal - cannot continue"<<endl;
                exit(-1);
            }
            else
            {
                set_ensemble_scale(pmdv,newscale);
                WriteTimeWindow(pmdv,tw_relative,outall);
                outall.close();
            }
            show_scaling_parameters(pmdv);
            if(engine_mode)
                ques="y";
            else
            {
                cout << "Answer y to accept this scale "
                    << "(anything else to try again): "
                    <<endl;
                cin >> ques;
                if(ques!="y")
                {
                    cout << "Enter new scale factor->";
                    cin >> newscale;
                }
            }
        }while(ques!="y");
        /* Now we write the string of files for animation.   
           This procedure does that */
        write_animation_files(pmdv,animate_twlength,animate_dt,
                animate_file_directory,animate_file_base,pvd_file,
                engine_mode);
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
