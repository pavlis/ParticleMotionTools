#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include "PMTimeSeries.h"
#include "seispp.h"
#include "dbpp.h"
#include "ThreeComponentSeismogram.h"
#include "PfStyleMetadata.h"
using namespace SEISPP;
void save_pmts(PMTimeSeries& d,string dir, string dfile_base,int band)
{
    const string base_error("Error in save_pmts procedure:  ");
    try {
        string full_fname;
        string sta=d.get_string("sta");
        long int evid=d.get_long("evid");
        stringstream ss;
        ss << dir <<"/"<<dfile_base<<"_"<<sta<<"_"<<evid<<".pmts";
        ofstream ofp;
        full_fname=ss.str();
        ofp.open(full_fname.c_str(),ios::out);
        if(ofp.fail()) throw SeisppError(base_error
                + "open failed for ofstream = "+full_fname);
        boost::archive::text_oarchive oa(ofp);
        oa << d;
        ofp.close();
    }catch(...){throw;};
}
bool dt_ok(ThreeComponentSeismogram& d,double target_dt,double tolerance)
{
    double ddt=fabs(d.dt-target_dt);
    ddt /= target_dt;
    if(ddt<tolerance)
        return true;
    else
        return false;
}
void build_working_view(DatascopeHandle& dbh, string subset_string)
{
    try{
        cout << "Building database working view"<<endl;
        DatascopeHandle ljhandle(dbh);
        dbh.lookup("event");
        dbh.natural_join("origin");
        dbh.subset("orid==prefor");
        dbh.natural_join("assoc");
        dbh.natural_join("arrival");
        cout << "Catalog view table size="<<dbh.number_tuples()<<endl;
        ljhandle.lookup("wfprocess");
        ljhandle.natural_join("evlink");
        ljhandle.natural_join("sclink");
        cout << "Number of waveform entries in wfprocess="
            << ljhandle.number_tuples()<<endl;
        list<string> jk;
        jk.push_back("evid");
        jk.push_back("sta");
        dbh.join(ljhandle,jk,jk);
        jk.clear();
        jk.push_back("sta");
        dbh.join(string("site"),jk,jk);
        cout << "Working table size="<<dbh.number_tuples()<<endl;
        if(subset_string.size()>0 && subset_string!="none")
        {
            dbh.subset(subset_string);
            cout << "Working table size after applying condition "
                << subset_string << "="<<dbh.number_tuples()<<endl;
        }
    }catch(...){throw;};
}
void usage()
{
    cerr << "mwpm db [-s subset -pf pffile]"<<endl
        << "db is assumed to contain 3c data output of extract_event program"
        <<endl
        << "Use -s to subset working view"<<endl;;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    /* As usual first parse the command line */
    if(argc<2) usage();
    int i,j;
    string dbname(argv[1]);
    string sstr("none");
    string pffile("mwpm.pf");
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-s")
        {
            ++i;
            if(i>=argc) usage();
            sstr=string(argv[i]);
        }
        else if (sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pffile=string(argv[i]);
        }
        else
            usage();
    }
    try {
        PfStyleMetadata control=pfread(pffile);
        /* Extract control parameters for this run first because 
           missing data and inconsistencies are the most common
           user error */
        MetadataList datamdl=get_mdlist(control,"input_data_mdl");
        /* Have to be dogmatic in this program that all data have
           a common sample rate.   This pair determines required
           sample rate and a tolerance factor (330s for example 
           will skew sample rate) */
        double target_dt=control.get_double("target_dt");
        double dt_tolerance=control.get_double("dt_fraction_error_tolerance");
        bool align_with_t0=control.get_bool("align_with_t0");
        string alignment_phase;
        if(align_with_t0)
            alignment_phase="t0";
        else
            alignment_phase=control.get_string("alignment_phase");
        double t0_offset(0.0);
        if(align_with_t0)
        {
            t0_offset=control.get_double("t0_offset");
            if(t0_offset<0)
            {
                cerr << "Illegal input t0_offset="<<t0_offset<<endl
                    << "Parameter must be positive"<<endl;
                exit(-1);
            }
        }
        /* This may change, but for now we always use this for 
           the keyword in metadata for defining t0*/
        const string alignkey("arrival.time");  
        /* This defines the window around t0 or alignment_phase
           that will be processed for particle motion estimates*/
        double cts,cte;
        cts=control.get_double("cut_time_window_start");
        cte=control.get_double("cut_time_window_end");
        TimeWindow cutwindow(cts,cte);
        /* Now define the multiwavelet transform operator object.
         Here we reuse the input parameter file.  The pf will contain
         a mix of data that way, but keeps the full processing 
         parameter set together.*/
        MWTransform mwt(pffile);
        int nbands=mwt.number_frequencies();
        /* These parameters define how and where serialized files
           are stored.
         Note file names will be constucted from dir+obname+ key
           string build from evid and sta names (evid and sta will thus 
           be required metadata */
        string outdir=control.get_string("output_data_directory");
        string obname=control.get_string("output_file_base_name");
        /* Particle motion ellipse algorithm has two fundamentall
           different approaches.   sample by sample method is assume
           if avlen is 1 or less.   */
        int avlen=control.get_int("particle_motion_time_average_length");
        int pmdt(1);
        if(avlen>1)
            pmdt=control.get_int("particle_motion_sampling_decimation_factor");

        AttributeMap am("css3.0");
        DatascopeHandle dbh(dbname,true);
        /* This acts like a subroutine and alters dbh to 
           become the handle to the working view */
        build_working_view(dbh,sstr);
        cout << "mwpm:  processing begins on database "<<dbname;
        if(sstr=="none")
            cout <<endl;
        else
            cout <<"using subset condition="<<sstr<<endl;
        long nrows=dbh.number_tuples();
        cout << "Number of rows in input view="<<nrows<<endl;

        for(dbh.rewind(),i=0;i<nrows;++i,++dbh)
        {
            try{
            auto_ptr<ThreeComponentSeismogram>
                d(new ThreeComponentSeismogram(dbh,datamdl,am));
            /* This helper (above) returns false if the data
               sample rate is out of tolerance */
            if(dt_ok(*d,target_dt,dt_tolerance))
            {
                if(align_with_t0)
                d->put(alignkey,d->t0+t0_offset);
                /* We always use this procedure that converst the 
                data to relative time using a phase as a reference.
                Note the logic above allows this to be a reference
                time defined by phase=T0.*/
                d=ArrivalTimeReference(*d,alignkey,cutwindow);
                MWTBundle dtransformed(*d,mwt);
                for(j=0;j<nbands;++j)
                {
                    PMTimeSeries pmts;
                    if(avlen>1)
                        pmts=PMTimeSeries(dtransformed,j,pmdt,avlen);
                    else
                        pmts=PMTimeSeries(dtransformed,j);
                    save_pmts(pmts,outdir,obname,j);
                }
            }
            else
            {
                cerr << "WARNING (mwpm):   data for station "
                    << d->get_string("sta")
                    << " has sample interval="<<d->dt 
                    <<" inconsistent with required dt="<<target_dt
                    <<endl
                    <<"Data from this station skipped"<<endl;
            }
            }catch(SeisppError& serr)
            {
                cerr << "Warning(mwpm):  data skipped for tuple "
                    << i << " of input view"<<endl
                    <<" when the following error was thrown"<<endl;
                serr.log_error();
            }
       }
    }catch(SeisppError& serr)
    {
        serr.log_error();
        exit(-1);
    }
    catch(std::exception& ex)
    {
        cerr << ex.what()<<endl;
        exit(-1);
    }
    catch(...)
    {
        cerr << "Something threw an unresolved exception"<<endl;
    }
}
