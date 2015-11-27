#include <stdio.h>
/*#include "su.h"*/
#include "par.h"
/* This is from su.h.  su.h generated a lot of errors so I pulled this
   in selectively*/
typedef union { /* storage for arbitrary type */
	char s[8];
	short h;
	unsigned short u;
	long l;
	unsigned long v;
	int i;
	unsigned int p;
	float f;
	double d;
	unsigned int U:16;
	unsigned int P:32;
} Value;
/* From newer version of su.h */
typedef char *cwp_String;
#include "segy.h"
#include "TimeSeries.h"
#include "HeaderMap.h"
#include "AttributeCrossReference.h"
TimeSeries ReadSegyTrace(FILE *fp, HeaderMap& hm, 
        AttributeCrossReference& xref, MetadataList& mdl)
{
    const string base_error("ReadSegyTrace:  ");
    segy tr;
    int iret;
    /*This is the base SU routine to read an SU trace C struct.
      It returns a nonzero value when successful and 0 if the read
      fails.  This is a bit obnoxious since it does not distinguish
      between a fail and an eof.   We handle this in a less than
      perfect way be returning an empty TimeSeries object if the 
      read failed. */
    iret=fgettr(fp,&tr);
    if(iret==0) return TimeSeries();
    TimeSeries d(tr.ns);
    /* fgettr converts xdr data to the segy interal struct.  Now 
       we use the HeaderMap to pull these out.   Note a less general
       solution worth considering is to pull out a frozen set of 
       names and convert them here to mesh with other particle motion
       code.   Will do it the more general way as this routine could
       have other uses. The code here is adapted from GenericFileHandle
     method LoadMetadata and the get template*/
    try {
        string extkey;
        MDtype keydatatype;
        MetadataList::iterator mdlptr;
        for(mdlptr=mdl.begin();mdlptr!=mdl.end();++mdlptr)
        {
            extkey=xref.external(mdlptr->tag);
            keydatatype=mdlptr->mdt;
            AttributeType rdtype=hm.dtype(extkey);
            short sival;
            int ival;
            double dval;
            string sval;
            bool bval;
            switch(keydatatype)
            {
                /* A header by definition starts at byte 0 so
                   we use the fact that the tr struct returned by
                   fgettr is a binary blob with a header. Passed
                   to HeaderMap method here as an opaque pointer*/
                case MDint:
                    if(rdtype==INT16)
                        sival=hm.get<short>(extkey,reinterpret_cast<unsigned char*>(&tr));
                    else
                        ival=hm.get<int>(extkey,reinterpret_cast<unsigned char*>(&tr));
                    d.put(mdlptr->tag,ival);
                    break;
                case MDreal:
                    dval=hm.get<double>(extkey,reinterpret_cast<unsigned char*>(&tr));
                    d.put(mdlptr->tag,dval);
                    break;
                case MDboolean:
                    bval=hm.get_bool(extkey,reinterpret_cast<unsigned char*>(&tr));
                    d.put(mdlptr->tag,bval);
                    break;
                case MDstring:
                    sval=hm.get_string(extkey,reinterpret_cast<unsigned char*>(&tr));
                    d.put(mdlptr->tag,sval);
                    break;
                case MDinvalid:
                default:
                    //Silently do nothing if invalid or something else
                    continue;
            }
        }
        /* We hard code these required attributes */
        d.ns=tr.ns;
        d.dt = ((double)tr.dt)*1.0e-6;
        if(tr.trid==2)
            d.live=false;
        else
            d.live=true;
    }catch(SeisppError& serr)
    {
        cerr << base_error <<"Error parsing header data."<<endl
            << "Header parsing loop threw this message:"<<endl;
        serr.log_error();
        cerr<<"Fatal Error: aborting"<<endl;
        exit(-1);
    }
    /* Now load the sample data - requires a float to double
       conversion */
    int i;
    for(i=0;i<tr.ns;++i) d.s.push_back((double)tr.data[i]);
    return d;
}


