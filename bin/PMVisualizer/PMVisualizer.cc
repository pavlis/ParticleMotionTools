#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <sstream>
#include "stock.h"
#include "seispp.h"
#include "Metadata.h"
#include "dbpp.h"
#include "PMVisualizerGUI.h"
using namespace std;
using namespace SEISPP;
string antelope_contrib_root_pf()
{
    string result;
    char *ant=getenv("ANTELOPE");
    if(ant==NULL) 
        result=string("./");
    else
    {
        result=string(ant);
        free(ant);
        result=result+"/contrib/data/pf";
    }
    return result;
}
/* This is the read routine for event data */
TimeSeriesEnsemble load_event_data(DatascopeHandle& dbh, string sstr,
        long evid, MetadataList& emdl, MetadataList& tmdl, AttributeMap& am)
{
    try{
        /* Dogmatically use wfprocess for this read approach 
         to match ParticleMotionVTKConverter comparable read routine*/
        dbh.lookup("evlink");
        stringstream ss;
        ss << "evid == "<<evid;
        dbh.subset(ss.str());
        /* Second subset passed */
        dbh.subset(sstr);
        dbh.natural_join("wfprocess");
        dbh.natural_join("event");
        dbh.natural_join("origin");
        dbh.subset(string("orid==prefor"));
        dbh.natural_join("assoc");
        dbh.natural_join("arrival");
        dbh.natural_join("site");
        int nmembers=dbh.number_tuples();
        if(nmembers<=0) throw SeisppError(string("load_event_data procedure:  ")
                + "No data found for requested event with subset condition "
                + sstr);
        TimeSeriesEnsemble d;
        d.member.reserve(3*nmembers);
        /* this is ugly - copied from database constructor.
        Point is to create ensemble metadata */
        Metadata ensmd(dbh,emdl,am);
        copy_selected_metadata(ensmd,dynamic_cast<Metadata&>(d),emdl);
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
            /* Here we need to take the 3C bundles apart.   Normally
               this function should only return 3 seismograms for
               the 3 channels */
            TimeSeries *chan;
            for(int k=0;k<3;++k)
            {
                chan=ExtractComponent(*s,k);
                d.member.push_back(*chan);
                delete chan;
            }
            delete s;
        }
        return d;
    }catch(...){throw;};
}
/* Time window based data loader.   Intentionally renamed the last
   arg to subset_channels instead of ensemble_plot_mode as the 
   late name is confusing in this context.   When true the subset
   string is passed as the chan subset condition.  Otherwise passed
   as sta condition */
TimeSeriesEnsemble load_time_window(DatascopeHandle& dbh, string sstr,
        TimeWindow tw,bool subset_channels)
{
    try {
        TimeSeriesEnsemble d;
        if(subset_channels)
            d=TimeSeriesEnsemble(dbh,tw,string("none"),sstr,true,true,false);
        else
            d=TimeSeriesEnsemble(dbh,tw,sstr,string("none"),true,true,false);
        return d;
    }catch(...){throw;};
}

/* This routine execs the ParticleMtionVTKConverter program in engine mode
   assembling all the pieces needed to make it run consistently with the gui.

Arguments:
    dbname - Datascope db name passed to ParticleMotionVTKConverter
    tw - time window (relative if event_mode is true, absolute otherwise)
    vtkparams - parameters returned by gui interface
    event_mode - critical switch controlling behaviour.  when true use  times
       defined by a phase arrival for event id evid
    evid - event id in event mode (ignored if event_mode is false)

Throws a SeisppError exception for several possible common failures.
*/
void run_vtk_converter(string dbname, TimeWindow tw, Metadata vtkparams,
        bool event_mode,long evid)
{
    const string base_error("run_vtk_converter procedure:  ");
    const string output_pffile("ParticleMotionVTKConverter.pf");
    try{
        /* Window may be absolute or relative, but for both modes we 
           set the internal parameter used by ParticleMotionVTKConverter even
           though a present this ony matters for event_mode.   Useful to 
           preserve time in either case in pf file used as input to that program*/
        vtkparams.put("process_window_start_time",tw.start);
        vtkparams.put("process_window_end_time",tw.end);
        /* Write results to a pf in the local directory.  This approach assumes
           ParticleMotionVTKConverter uses the antelope pf search path to find
           default parameters that are not set by the gui. */
        ofstream opf;
        opf.open(output_pffile.c_str());
        if(opf.fail())
            throw SeisppError(base_error
                    + "open failed on Pf file="+output_pffile);
        opf << vtkparams;
        opf.close();
        /* We will use execlp which we assume will find ParticleMotionVTKConverter
           in the shell path. */
        const string progname("ParticleMotionVTKConverter");
        int iret;  // used to hold execlp return code 
        if(event_mode)
        {
            stringstream ss;
            ss << evid;
            string evstr=ss.str();
            iret=execlp(progname.c_str(),progname.c_str(),dbname.c_str(),
                    "-e",evstr.c_str(),
                    "-engine",NULL);
        }
        else
        {
            iret=execlp(progname.c_str(),progname.c_str(),dbname.c_str(),
                    "-a",strtime(tw.start),strtime(tw.end),
                    "-engine",NULL);
        }
        if(iret) throw SeisppError(base_error
                + "execlp failed running ParticleMotionVTKConverter");
    }catch(...){throw;};
}

void usage()
{
    cerr << "PMVisualizer db (-e evid | -ts t0 ) (-3C sta | -comp XXX) [-pf pffile]" <<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    /* First as usual crack the command line. */
    if(argc<6) usage();
    string dbname(argv[1]);
    long evid(-1);
    double start_time(0.0);
    bool event_mode(true);
    bool ensemble_plot_mode(true);
    bool data_mode_set(false);
    bool plot_mode_set(false);
    string comp,sta;
    string sarg;
    int i;
    for(i=2;i<argc;++i)
    {
        sarg=string(argv[i]);
        if(sarg=="-comp")
        {
            ensemble_plot_mode=true;
            ++i;
            if(i>=argc) usage();
            comp=string(argv[i]);
            plot_mode_set=true;
        }
        else if(sarg=="-3C")
        {
            ensemble_plot_mode=false;
            ++i;
            if(i>=argc) usage();
            sta=string(argv[i]);
            plot_mode_set=true;
        }
        else if(sarg=="-e")
        {
            event_mode=true;
            ++i;
            if(i>=argc) usage();
            evid=atol(argv[i]);
            data_mode_set=true;
        }
        else if(sarg=="-ts")
        {
            event_mode=false;
            ++i;
            if(i>=argc) usage();
            start_time=str2epoch(argv[i]);
            data_mode_set=true;
        }
        else
            usage();
    }
    if(!data_mode_set)usage();
    if(!plot_mode_set)usage();
    string cr=antelope_contrib_root_pf();
    string pffile=cr+"PMVisualizer.pf";
    for(i=5;i<argc;++i)
    {
        sarg=string(argv[i]);
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
        Pf *pf;
        if(pfread(const_cast<char *>(pffile.c_str()),&pf))
        {
            cerr << argv[0]<<": pfread failed on file="<<pffile<<endl;
            exit(-1);
        }
        Metadata plotmd(pf);
        /* These can and should be identical to ParticleMotionVTKConverter*/
        MetadataList tracemdl=pfget_mdlist(pf,"trace_mdlist");
        MetadataList ensemblemdl=pfget_mdlist(pf,"ensemble_mdlist");
        AttributeMap am("css3.0");
        DatascopeHandle dbh(dbname,true);
        string subset_str;
        if(ensemble_plot_mode)
            subset_str="chan=~/"+comp+"/";
        else
            subset_str="sta=~/"+sta+"/";

        /* Different data options are implemented as different procedures
           that all return a TimeSeriesEnsemble.  Done
           to simplify logic of main.   Use a copy instead of a pointer
           as the size of d should not be huge */
        TimeSeriesEnsemble d;
        if(event_mode)
            d=load_event_data(dbh,subset_str,evid,ensemblemdl,tracemdl,am);
        else
        {
            double endtime=start_time+plotmd.get_double("DisplayWindowLength");
            d=load_time_window(dbh,subset_str,TimeWindow(start_time,endtime),
                        ensemble_plot_mode);
        }
        PMVisualizerGUI  win(plotmd);
        char ques;
        do{
            do{
                /* This pf is the output of the gui */
                const string gui_pffile("PMVisualizerGUI.pf");
                cout << "Fill out any changes in the tk GUI window and push accept"
                    <<endl;
                /* probably need a fifo to receive output of accept button */
                PfStyleMetadata guipf=pfread(gui_pffile);
                win.load_pfdata(guipf);
                win.filter_and_plot(d);
                cout << "Enter y to use this filter.   "
                    << "Type any other key to try again:";
                ques=getchar();
            }while(ques=='y');
            cout << "Select the time window for display in paraview"<<endl;
            TimeWindow tw=win.get();
            Metadata vtkparams=win.get_parameters();
            run_vtk_converter(dbname,tw,vtkparams,event_mode,evid);
            cout << "Check output in paraview"<<endl;
            cout << "Try again?"<<endl
                << "If so enter y, type any other key to exit"<<endl;
            ques=getchar();
        }while(ques=='y');

    }catch(SeisppError& serr)
    {
        serr.log_error();
        exit(-1);
    }
    catch(std::exception& stdexcp)
    {
        cerr << stdexcp.what();
        exit(-1);
    }
}





