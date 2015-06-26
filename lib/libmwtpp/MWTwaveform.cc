#include "MWTransform.h"
MWTwaveform::MWTwaveform(MWtrace& mw) : ComplexTimeSeries(mw.nz)
{
    /* first set attributes of an MWTwaveform */
    dt0=mw.dt0;
    decimation_factor=mw.decimation_factor;
    f0=mw.f0;
    fw=mw.fw;
    /* Now set attributes in mw shared by ComplexTimeSeries*/
    ns=mw.nz;
    dt=mw.dt;
    t0=mw.starttime;
    // Use of epoch time is implicity in MWtrace we force it 
    tref=absolute;
    /* note we drop endtime - a method in BasicTimeSeries */
    int i;
    /* MWtrace uses 32 bit float complex values.  We have to 
       take them apart and put them into C++ form */
    for(i=0;i<ns;++i)
    {
        float re=mw.z[i].r;
        float im=mw.z[i].i;
        Complex zcomp((double)re,(double)im);
	/* We use the indexing operator instead of push because
	the ComplexTimeSeries constructor initializes an initial
	vector of length ns to zeros */
	s[i]=zcomp;
    }
    live=true;
    /* finally we push the private attributes to metadata to 
       make them accessible.  Avoids excessive getters */
    put("dt0",dt0);
    put("decimation_factor",decimation_factor);
    put("f0",f0);
    put("fw",fw);
}
MWTwaveform::MWTwaveform(const MWTwaveform& parent) : ComplexTimeSeries(parent)
{
    dt0=parent.dt0;
    decimation_factor=parent.decimation_factor;
    f0=parent.f0;
    fw=parent.fw;
}
MWTwaveform& MWTwaveform::operator=(const MWTwaveform& parent)
{
    if(this!=&parent)
    {
        this->ComplexTimeSeries::operator=(parent);
        dt0=parent.dt0;
        decimation_factor=parent.decimation_factor;
        f0=parent.f0;
        fw=parent.fw;
    }
    return(*this);
}
