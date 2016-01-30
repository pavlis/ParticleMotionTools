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
class MDTable
{
    public:
        MDTable(PfStyleMetadata pfmd);
        int number_attributes()
        {
            return(nm.size());
        };
        int number_tuples()
        {
            return(strval.size());
        };
        /*! Cautiously set a field from the table.

          This is the primary working method of this object.
          Intended to be applied using a loop over attributes
          without needing the baggage of type cracking. Will
          set attribute at position j associated with tuple i. 
          
          \param d - Metadata object to set
          \param  i - tuple to access
          \param  j - set value associated with attribute at column j

          \return - name of attribute set
          \exception - throws a SeisppError if i and j are outside
          bounds of the table
          */
        string set(Metadata& d,int i, int j);
        /*! Set all attributes defined in tuple j.
        
        This method should be used if an entire list of attributes
        is always set.  
        
        \param d - Metadata object to set
        \param i - tuple of attributes to be set. 

        \return number of attributes set
        */
        int set(Metadata& d,int i);
        /*! Return name of attribute at position j */
        string name(int j){return nm[j];};
        /*! Return type of attribute at position j */
        MDtype type(int j){return t[j];};
    private:
        vector<string> nm;
        vector<MDtype> t;
        vector<string> strval;
};
MDTable::MDTable(PfStyleMetadata pfmd)
{
    MetadataList mdl=get_mdlist(pfmd,"types");
    MetadataList::iterator mptr;
    for(mptr=mdl.begin();mptr!=mdl.end();++mptr)
    {
        this->nm.push_back(mptr->tag);
        this->t.push_back(mptr->mdt);
    }
    list<string> vlist=pfmd.get_tbl("values");
    /* interface returns data we need as a list.  Need to 
       convert to a vector container*/
    list<string>::iterator vptr;
    for(vptr=vlist.begin();vptr!=vlist.end();++vptr)
        strval.push_back(*vptr);
}
string MDTable::set(Metadata& d,int i0, int j0)
{
    const string range_error("MDTable::set method:  index out of range\n");
    if(i0<0) throw SeisppError(range_error + "tuple index requested was negative");
    if(j0<0) throw SeisppError(range_error 
            + "attribute (column) index requested was negative");
    if(i0>=strval.size()) throw SeisppError(range_error
            + "tuple index requested is larger than table size");
    if(j0>=nm.size()) throw SeisppError(range_error
            + "attribute (column) index requested is larger than list of attributes");
    stringstream ss(strval[i0]);
    int j;
    for(j=0;j<=j0;++j)
    {
        double dval;
        long ival;
        string s_now;
        ss >> s_now;
        if(j<j0) continue;
        switch(this->t[j])
        {
            case MDreal:
                dval=atof(s_now.c_str());
                d.put(this->nm[j],dval);
                break;
            case MDint:
                ival=atol(s_now.c_str());
                d.put(this->nm[j],ival);
                break;
            case MDstring:
                d.put(this->nm[j],s_now);
                break;
            case MDinvalid:
            default:
                throw SeisppError(string("MDtable::set method")
                        + " Parsing problem.  Name="
                        + this->nm[j] + " has invalid type defined");
        }
    }
    return(nm[j]);
}


int MDTable::set(Metadata& d,int i0)
{
    const string range_error("MDTable::set method:  index out of range\n");
    if(i0<0) throw SeisppError(range_error + "tuple index requested was negative");
    if(i0>=strval.size()) throw SeisppError(range_error
            + "tuple index requested is larger than table size");
    stringstream ss(strval[i0]);
    int j,count;
    for(j=0,count=0;j<this->number_attributes();++j)
    {
        double dval;
        long ival;
        string str;
        switch(this->t[j])
        {
            case MDreal:
                ss>>dval;
                d.put(this->nm[j],dval);
                ++count;
                break;
            case MDint:
                ss>>ival;
                d.put(this->nm[j],ival);
                ++count;
                break;
            case MDstring:
                ss>>str;
                d.put(this->nm[j],str);
                ++count;
                break;
            case MDinvalid:
            default:
                cerr << "MDtable::set(WARNING):  "
                    << "Attribute "<<this->nm[j]
                    << " has invalid type.  Attribure NOT SET"<<endl;
        }
    }
    return count;
}
void usage()
{
    cerr << "set_metadata [-pf pffile] < infile > outfile"
        <<endl
        << "Input and output are serialized ThreeComponentEnsemble objects"
        <<endl
        << "sets ensemble and/or member metadata using a pf format"<<endl
        << "Assumes members of ensemble are 3C seismogram objects"
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
    string pffile("set_metadata.pf");
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
        const string emdkey("EnsembleMetadata");
        const string memkey("MemberMetadata");
        PfStyleMetadata control=pfread(pffile);
        /* The concept here is that this object will contain
           metadata to be set as global for the ensemble */
        PfStyleMetadata ensmd=control.get_branch(emdkey);
        MDTable ensmdtbl(ensmd);
        /* These entries will have a vector of data for the ensemble 
           with a tag */
        PfStyleMetadata memmd=control.get_branch(memkey);
        MDTable memmdtbl(memmd);
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentEnsemble d;
        d=read_object<ThreeComponentEnsemble>(ia);
        /* This is made a fatal error for now.  It perhaps should
           be handled more gracefully.  Certainly would need to be if
           this were in a larger system */
        if(memmdtbl.number_tuples() != d.member.size())
        {
            cerr << "set_metadata(FATAL ERROR):  "
                << "Size mistmatch in table in pffile="<<pffile
                << " and input data file"<<endl
                << "Table size="<<memmdtbl.number_tuples()
                << " while number of seismograms in ensemble = "
                << d.member.size()<<endl;
            exit(-1);
        }
        /* Sets this ensemble metadata  */
        ensmdtbl.set(dynamic_cast<Metadata&>(d),0);
        int count_returned,count_expected;
        count_expected=memmdtbl.number_attributes();
        for(i=0;i<d.member.size();++i)
        {
            count_returned=memmdtbl.set(dynamic_cast<Metadata&>(d.member[i]),i);
            if(count_returned!=count_expected)
                cerr << "Warning - some attributes not set for member "
                    <<i<<endl;
        }
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

