#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <list>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "seispp.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
/* You will get lots of errors without these namespace
   declaration*/
using namespace std;   // most compilers do not require this
using namespace SEISPP;  //This is essential to use SEISPP library
/* You shoudl always include a usage procedure like this to trap
   command line parsing problems. */
void usage()
{
    cerr << "gather -i key1 key2 ... -s key1 key2 ... "
        <<endl
        << "Build gathers using list of integer and string keys"<<endl
        << "Follow -i with list of integer keys that define a match"<<endl
        << "Follow -s with a list of string keys that define a match"<<endl
        << "Gather is a group by requiring exactly matching all key values"<<endl
        << "group keys are posted to ensemble metadata"<<endl;
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
bool keys_match(ThreeComponentSeismogram& d,list<string>& sk, list<string>& ik,
    list<string>stest, list<int> itest)
{
  try{
    list<string>::iterator skptr,svptr;
    for(skptr=sk.begin(),svptr=stest.begin();skptr!=sk.end();++skptr,++svptr)
    {
      string svnow=d.get_string(*skptr);
      if(svnow != (*svptr)) return false;
    }
    list<string>::iterator ikptr;
    list<int>::iterator ivptr;
    for(ikptr=ik.begin(),ivptr=itest.begin();ikptr!=ik.end();++ikptr,++ivptr)
    {
      int intvnow=d.get_int(*ikptr);
      if(intvnow != (*ivptr)) return false;
    }
    return true;
  }catch(...){throw;};
}
void reset_ensemble(ThreeComponentSeismogram& d,list<string> skeys,
  list<string> ikeys,ThreeComponentEnsemble& dout,
    list<string>& svaltest,list<int>& ivaltest)
{
  try{
    list<string>::iterator sptr;
    list<string>::iterator iptr;
    dout.member.clear();
    svaltest.clear();
    ivaltest.clear();
    for(sptr=skeys.begin();sptr!=skeys.end();++sptr)
    {
      string sval=d.get_string(*sptr);
      svaltest.push_back(sval);
    }
    for(iptr=ikeys.begin();iptr!=ikeys.end();++iptr)
    {
      int ival=d.get_int(*iptr);
      ivaltest.push_back(ival);
    }
  }catch(...){throw;};
}

bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    int i;
    if(argc<3)usage();
    list<string> ikeys,skeys;
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-i")
        {
            ++i;
            if(i>=argc)break;
            do
            {
              sarg=argv[i];
              ikeys.push_back(sarg);
              ++i;
            }while(i<argc && (sarg!="-s"));
            if(sarg=="-s")
            {
              --i;
              continue;
            }
        }
        else if(sarg=="-s")
        {
          ++i;
          if(i>=argc)break;
          do
          {
            sarg=argv[i];
            skeys.push_back(sarg);
            ++i;
          }while(i<argc && (sarg!="-i"));
          if(sarg=="-i")
          {
            --i;
            continue;
          }
        }
        else
            usage();
    }
    try{
        list<string> svaltest;
        list<int> ivaltest;
        boost::archive::text_iarchive ia(cin);
        boost::archive::text_oarchive oa(cout);
        ThreeComponentSeismogram d;
        ThreeComponentEnsemble dout;
        list<string>::iterator sptr,iptr;
        int nseis,ngather;
        nseis=0;   ngather=0;
        do{
          try{
            ia >> d;
            if(nseis==0)
            {
              for(sptr=skeys.begin();sptr!=skeys.end();++sptr)
              {
                string stmp=d.get_string(*sptr);
                svaltest.push_back(stmp);
                dout.put(*sptr,stmp);
              }
              for(iptr=ikeys.begin();iptr!=ikeys.end();++iptr)
              {
                int itmp=d.get_int(*iptr);
                ivaltest.push_back(itmp);
                dout.put(*iptr,itmp);
              }
            }
            else
            {
              if(keys_match(d,skeys,ikeys,svaltest,ivaltest))
              {
                dout.member.push_back(d);
              }
              else
              {
                ++ngather;
                oa << dout;
                /* This clears dout and then initializes ensemble metadata
                and sets svaltest and ivaltest with values from d.*/
                reset_ensemble(d,skeys,ikeys,dout, svaltest,ivaltest);
              }
            }
          }catch(boost::archive::archive_exception const& e)
          {
            ++ngather;
            oa << dout;
            /* Oddity of boost serialization implementation is this seems
            the only way to string objects together and have a data driven
            input.  Normal exit as an exception is evil, but there seems
            no choice. */
            cerr << "Gather processed "<<nseis<<" 3c seismograms"<<endl
                << "Assembled "<<ngather<<" ensembles"<<endl;
            exit(0);
          }
          ++nseis;
        }while(1);
    }catch(SeisppError& serr)
    {
        serr.log_error();
        exit(-1);
    }
    catch(std::exception& stexc)
    {
        cerr << stexc.what()<<endl;
        exit(-1);
    }
}
