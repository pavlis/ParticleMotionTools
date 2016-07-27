#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "seispp.h"
#include "ensemble.h"
#include "PMTimeSeries.h"
using namespace std;   // most compilers do not require this
using namespace SEISPP;  //This is essential to use SEISPP library
void usage()
{
    cerr << "listhdr infile [-csv format_file -t objecttype]  >outfile"
        <<endl
        << "List metadata components of a stream of serialized objects"
        <<endl
        << " infile is concatenated set of serialized objects"<<endl
        << " -csv - write the output as a csv file using the format defined in format_file"
        <<endl
        << "        (default dumps all with operator <<"<<endl
        << " -t - specify the type of object expected"<<endl
        << "      (Currently accept:  ThreeComponentSeismogram (default), ThreeComponentEnsemble, "<<endl
        << "and PMTimeSeries)"
        <<endl;;
    exit(-1);
}
enum AllowedObjects {TCS, TCE, PMTS};
/* Special read function for this program.  Will return metadata portion of any
   object that is a child of Metadat */
template <class InputObject> Metadata
    read_object(boost::archive::text_iarchive& ia)
{
    InputObject d;
    try{
        ia>>d;
    }catch(...)
    {
        throw SeisppError(string("read_object failed:  ")
                    + "Check that input file is a boost text archive file");
    }
    return Metadata(d);
}
template <class InputObject> vector<size_t> 
        build_index(ifstream& ifs)
{
    /* Safety valve to avoid a runaway */
    int nd(0);
    const int ndmax(1000000);
    vector<size_t> foff_values;
    boost::archive::text_iarchive ia(ifs);
    InputObject d;
    do{
      try{
        ia >> d;
        size_t foff=ifs.tellg();
        foff_values.push_back(foff);
        ++nd;
      }catch(boost::archive::archive_exception const& e)
      {
          break;
      }
    }while(nd<ndmax);  // safer than infiite loop
    return foff_values; 
}
/*! Simple class to drive csv outputs in this program. */
class MetadataComponent
{
    public:
        string key; //metadata key to access 
        MDtype mdt;  //type enum (defined in Metadata.h)
        string undefined_value;  // This value is value to write to output if key not found
};
vector<MetadataComponent> parse_csv_format_file(string fname)
{
    ifstream icff;
    icff.open(fname.c_str(),ios::in);
    if(!icff) 
    {
        cerr << "parse_csv_format_file procedure:  open failed on file "
            << fname<<endl;
        exit(-1);
    }
    vector<MetadataComponent> result;
    char inpline[128];
    while(icff.getline(inpline,128))
    {
        stringstream ss(inpline);
        MetadataComponent mc;
        ss>>mc.key;   
        string mdttest;
        ss >> mdttest;
        if( (mdttest=="double") || (mdttest=="MDreal") || (mdttest=="real") || (mdttest=="float") )
            mc.mdt=MDreal;
        else if( (mdttest=="int") || (mdttest=="MDint") || (mdttest=="log") )
            mc.mdt=MDint;
        else if(mdttest=="boolean")
            mc.mdt=MDboolean;
        else if( (mdttest=="string") || (mdttest=="String") || (mdttest=="MDstring") )
            mc.mdt=MDstring;
        else
        {
            cerr << "parse_csv_format_file procedure: Unrecognized type ="<<mdttest
               <<" for entry with key="<<mc.key<<endl
              << "Fatal error - exiting"<<endl;
           exit(-1);
        }
        ss>>mc.undefined_value;
        result.push_back(mc);
    } 
    return result;
}
/* This procedure will write one line to csv file extracting from Metadata d using
   the mdl vector to define the output format and null defaults */
int WriteToCSVFile(Metadata& d,ostream& ofs,vector<MetadataComponent>& mdl)
{
    int nlast=mdl.size()-1;
    int i;
    int nsaved;
    double dval;
    int ival;
    bool bval;
    string sval;
    vector<MetadataComponent>::iterator mdptr;
    /* We do this to be sure epoch times get printed correctly in all cases.   Makes
       sense here as bloated files are unlikely to be an issue with metadata output*/
    ofs<<std::setprecision(13);
    for(i=0,mdptr=mdl.begin();mdptr!=mdl.end();++mdptr,++i)
    {
        try{
            switch (mdptr->mdt)
            {
                case MDreal:
                    dval=d.get_double(mdptr->key);
                    ofs << dval;
                    break;
                case MDint:
                    ival=d.get_int(mdptr->key);
                    ofs << ival;
                    break;
                case MDstring:
                    sval=d.get_string(mdptr->key);
                    ofs << sval;
                    break;
                case MDboolean:
                    bval=d.get_bool(mdptr->key);
                    ofs << bval;
                    break;
                default:
                    cerr << "Illegal mdtype specified"<<endl
                        << "This should not happen unless the program overwrites itself"<<endl
                        << "Fatal error:  exiting"<<endl;
                    exit(-1);
            }
            ++nsaved;
        }catch (SEISPP::MetadataGetError& err)
        {
            ofs << mdptr->undefined_value;
        }
        if(i<nlast) ofs<<",";
    }
    ofs << endl;   // endl does not work for minicsv
    return nsaved;
}

AllowedObjects get_object_type(string otype)
{
    if(otype=="ThreeComponentSeismogram")
        return TCS;
    else if(otype=="ThreeComponentEnsemble")
        return TCE;
    else if(otype=="PMTimeSeries")
        return PMTS;
    else
    {
        cerr << "Do not know how to handle object type="<<otype
            <<endl<< "Cannot continue"<<endl;
        exit(-1);
    }
}



bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    int i;
    if(argc<2) usage();
    string infile(argv[1]);
    bool csv_output(false);
    string otype("ThreeComponentSeismogram");
    string fname_csvo;
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-csv")
        {
            ++i;
            if(i>=argc)usage();
            csv_output=true;
            fname_csvo=string(argv[i]);
        }
        else if(sarg=="-t")
        {
            ++i;
            if(i>=argc)usage();
            otype=string(argv[i]);
        }
        else
            usage();
    }
    try{
        AllowedObjects dtype=get_object_type(otype);
        ifstream ifs;
        ifs.open(infile.c_str(),ios::in);
        if(!ifs)
        {
            cerr << "Open failed on input data file"<<infile<<endl;
            exit(-1);
        }
        vector<MetadataComponent> csv_format_info;
        if(csv_output)
            csv_format_info=parse_csv_format_file(fname_csvo);
        boost::archive::text_iarchive ia(ifs);
        vector<size_t> fofflist;
        switch (dtype)
        {
            case TCS:
                fofflist=build_index<ThreeComponentSeismogram>(ifs);
                break;
            case TCE:
                fofflist=build_index<ThreeComponentEnsemble>(ifs);
                break;
            case PMTS:
                fofflist=build_index<PMTimeSeries>(ifs);
                break;
            default:
                cerr << "Coding problem - dtype variable does not match enum"
                    <<endl
                    << "Fatal error - bug fix required. "<<endl;
                exit(-1);
        };
        int nd=fofflist.size();
        for(i=0;i<nd;++i)
        {
            size_t foff=fofflist[i];
            ifs.seekg(foff,ifs.beg);
            Metadata d;
            switch(dtype)
            {
                case TCS:
                    d=read_object<ThreeComponentSeismogram>(ia);
                    break;
                case TCE:
                    d=read_object<ThreeComponentEnsemble>(ia);
                    break;
                case PMTS:
                    d=read_object<PMTimeSeries>(ia);
                    break;
                default:
                    cerr << "Unrecognized data type:  This should not happen"
                        << "Program has probably overwritten itself from a bug"
                        <<endl
                        << "Cannot continue.  Exiting."<<endl;
                    exit(-1);
            };
            if(csv_output)
            {
                WriteToCSVFile(d,cout,csv_format_info);
            }
            else
                cout << d;

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

