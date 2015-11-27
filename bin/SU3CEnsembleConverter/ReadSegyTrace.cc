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
using namespace SEISPP;
TimeSeries ReadSegyTrace(FILE *fp)
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
    TimeSeries d;
    d.s.reserve(tr.ns);
    /* We extract only a fixed set of header values.   This is not
       general, but the first attempt at this using a HeaderMap 
       object was more hassle than it was worth.  Instead this
       just hard codes what is extracted.   Hack this if 
       the code needs to be adapted to another project.  
       
       First the required elements for a TimeSeries object. */
    try{
        d.ns=tr.ns;
        d.dt = ((double)tr.dt)*1.0e-6;
        if(tr.trid==2)
            d.live=false;
        else
            d.live=true;
        d.ns=tr.ns;
        d.dt = ((double)tr.dt)*1.0e-6;
        if(tr.trid==2)
            d.live=false;
        else
            d.live=true;
        /* Now we copy a few things to the Metadata header */
        d.put("tracl",tr.tracl);
        d.put("tracr",tr.tracr);
        d.put("fldr",tr.fldr);
        d.put("tracf",tr.tracf);
        // alias for downstream
        char sbuf[10];
        sprintf(sbuf,"%d",tr.tracf);
        d.put("chan",sbuf);
        d.put("ep",tr.ep);
        // another useful alias
        d.put("evid",tr.ep);
        d.put("nvs",tr.nvs);
        /* These are coordinates.  Here we store the 
           raw values as ints and a real value 
           that will be actually used downstream */
        double dcoord;
        double scale=(double)(tr.scalco);
        dcoord=(double)(tr.gx);
        dcoord/=scale;
        d.put("rx",dcoord); // note conversion of gx,gy to rx,ry
        dcoord=(double)(tr.gy);
        dcoord/=scale;
        d.put("ry",dcoord);
        dcoord=(double)(tr.sx);
        dcoord/=scale;
        d.put("sx",dcoord);
        dcoord=(double)(tr.sy);
        dcoord/=scale;
        d.put("sy",dcoord);
        d.put("relev",tr.gelev);
        d.put("selev",tr.selev);
        /* This is a peculiar Homestake data oddity.  offset
           is in m and does not include the scalco factor.  This 
           may be segy standard, but beware as I'm not sure. */
        d.put("offset",(double)tr.offset);
        /* Double these in metadata */
        d.put("nsamp",(int)tr.ns);
        d.put("int_dt",(int)tr.dt);
        d.put("dt",d.dt);
        d.put("samprate",1.0/(d.dt));
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
    for(i=0;i<tr.ns;++i) 
        d.s.push_back((double)tr.data[i]);
    return d;
}


