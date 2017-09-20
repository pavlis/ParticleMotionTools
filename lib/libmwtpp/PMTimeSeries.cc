#include <math.h>
#include "PMTimeSeries.h"
using namespace SEISPP;
void ComputePMStats(vector<ParticleMotionEllipse>& d,
                ParticleMotionEllipse& avg, ParticleMotionError& err);
/* These come from multiwavelet.h.   Extracted only these to avoid symbol
   collisions if the entire file is included (complex number mismatch
   between an older C code and C++).  MAINTENANCE ISSUE - WARNING. 
 Note that ComputePMStatis is effectively an interface routine
 between these structs and the new C++ version*/

typedef struct Particle_Motion_Ellipse_{
	double major[3];
	double minor[3];
	double rectilinearity;  
} Particle_Motion_Ellipse;

typedef struct Particle_Motion_Error_{
	double dtheta_major,dphi_major;  /* errors in spherical coordinate
					angles for major ellipse direction */
	double dtheta_minor,dphi_minor;  /* same for minor axis */
	double delta_rect;  /* error in rectilinearity */
	int ndgf_major, ndgf_minor, ndgf_rect;  /* Degrees of freedom */
} Particle_Motion_Error;
extern "C" {
void pmvector_average(Particle_Motion_Ellipse *, int ,
              Particle_Motion_Ellipse *, Particle_Motion_Error *);
}
	
/* Private method used by constructors once attributes are set to put
   private data into metadata object.   No need then for getters in
   the interface */
void PMTimeSeries::post_attributes_to_metadata()
{
    this->put("averaging_length",averaging_length);
    this->put("f0",f0);
    this->put("fw",fw);
    this->put("decfac",decfac);
    /* These should normally be set in a valid time series, but best
       to be sure */
    this->put("nsamp",ns);
    this->put("samprate",1.0/dt);
}
PMTimeSeries::PMTimeSeries() : Metadata(), BasicTimeSeries()
{
    averaging_length=0;
    f0=0.0;
    fw=0.0;
    decfac=0;
}

PMTimeSeries::PMTimeSeries(MWTBundle& d, int band, int timesteps, int avlen)
    : Metadata(dynamic_cast<Metadata&>(d))
{
    const string base_error("PMTimeSeries time averaging constructor:  ");
    /* Some sanity checks */
    if(d.number_members() != 3) 
    {
        stringstream ss;
        ss << "Bundle size passed = "<<d.number_members()<<endl
            << "This must be exactly 3 - likely coding error"<<endl;
        throw SeisppError(base_error + ss.str());
    }
    if(band>=d.number_bands()) 
    {
        stringstream ss;
        ss << "Illegal request for band="<<band<<endl
            << "band requested must be between 0 and "<<d.number_bands()-1<<endl;
        throw SeisppError(base_error+ss.str());
    }
    int i,iw,it;
    try {
        /* First set the private scalar attributes, then move to 
           build the large data vectors */
        averaging_length=avlen;
        f0=d.get_f0(band);
        fw=d.get_fw(band);
        //dt=d.sample_interval(band);
        decfac=d.get_decfac(band)*timesteps;
        int iw,i;
        double up[3]={0.0,0.0,1.0};
        vector<MWTwaveform> x,y,z;
        int nw=d.number_wavelets();
        x.reserve(nw);  y.reserve(nw);   z.reserve(nw);
        /* First load the component waveforms for this band in the x,y, and z vectors. */
        for(iw=0;iw<nw;++iw)
        {
            x.push_back(d(band,iw,0));
            y.push_back(d(band,iw,1));
            z.push_back(d(band,iw,2));
        }
        /* These are required attributes from BasicTimeSeries.  We have
           to assume we can estract them from the components we just
           built. */
        double dtparent;
        dtparent = x[0].get_dt0();
        this->dt=(double)dtparent*((double)decfac);;
        this->ns=(x[0].ns)/decfac;
        this->tref=x[0].tref;
        this->t0=x[0].t0;
        this->live=false;   // Set this way in case we throw an error
        /* These floats are in units of s - derived from steps in samples */
        double time_avlen=(dtparent)*((double)(avlen));
        /* This is a sanity check.   Perhaps unnecessary overhead, but
         cost is not high for stability it can buy */
        for(iw=0;iw<nw;++iw)
        {
            if((x[iw].ns!=y[iw].ns) || (y[iw].ns!=z[iw].ns))
            {
                stringstream ss;
                ss << "Size mismatch in ComplexTimeSeries components for wavelet="
                    << iw<<endl;
                throw SeisppError(base_error + ss.str());
            }
            if((x[iw].t0!=t0) || (y[iw].t0!=t0) || (z[iw].t0!=t0) )
            {
                stringstream ss;
                ss << "Start time mismatch in ComplexTimeSeries components for wavelet="
                    << iw<<endl;
                throw SeisppError(base_error + ss.str());
            }
        }
        /* This traps a condition that can easily happen at 
           low frequencies if the original time window is too short */
        int nsbandrequested;
        nsbandrequested=(int)((x[0].endtime()-x[0].t0)/x[0].dt);
        if(avlen>nsbandrequested) 
        {
            stringstream ss;
            ss << "Averaging length too long for input data segment"<<endl
                << "Data segement in band "<<band<<" is only "<<nsbandrequested<<" samples long"<<endl
                << "This less than requested averaging length of "<<avlen<<endl;
            throw SeisppError(base_error + ss.str());
        }
        vector<ParticleMotionEllipse> pmi;   // nw estimates for each i
        ParticleMotionEllipse avg;
        ParticleMotionError err;
        pmi.reserve(nw);
        //assume x,y, and z have common start times 
        this->t0=x[0].t0+time_avlen/2.0;  //use centered time as reference
        //this->dt=timesteps;
        //cout << "x[0].endtime() is " << x[0].endtime() << endl;
        double t;  // this is stare time of averaging window not center
        for(i=0,t=(this->t0);i<ns;++i,t+=(this->dt))
        {
            TimeWindow tw(t,t+time_avlen);
            //cout << tw.end << endl;
            if ( tw.end < x[0].endtime() )
            {
                for(iw=0;iw<nw;++iw)
                    pmi.push_back(ParticleMotionEllipse(x[iw],y[iw],z[iw],tw,up));
                ComputePMStats(pmi,avg,err);
                pmdata.push_back(avg);
                pmerr.push_back(err);
                pmi.clear();
            }
            else { break; }
        }
        /* Reset ns if necessary.   Do this silently unless ns is 0 or less */
        if(pmdata.size()<=0) throw SeisppError(base_error
                + "Data window is too short for specified parameters - zero length PMTimeSeries result");
        if((this->ns) != pmdata.size()) this->ns = pmdata.size();
        this->post_attributes_to_metadata();
        live=true;
    }catch(...){throw;};
}
PMTimeSeries::PMTimeSeries(MWTBundle& d, int band)
    : Metadata(dynamic_cast<Metadata&>(d))
{
    const string base_error("PMTimeSeries sample-by-sample constructor:  ");
    if(band>=d.number_bands()) 
    {
        stringstream ss;
        ss << "Illegal request for band="<<band<<endl
            << "band requested must be between 0 and "<<d.number_bands()-1<<endl;
        throw SeisppError(base_error+ss.str());
    }
    int i,iw,it;
    try {
        /* First set the private scalar attributes, then move to 
           build the large data vectors */
        averaging_length=1;
        f0=d.get_f0(band);
        fw=d.get_fw(band);
        dt=d.sample_interval(band);
        decfac=d.get_decfac(band);
        int iw,i;
        double up[3]={0.0,0.0,1.0};
        vector<MWTwaveform> x,y,z;
        int nw=d.number_wavelets();
        x.reserve(nw);  y.reserve(nw);   z.reserve(nw);
        /* First load the component waveforms for this band in the x,y, and z vectors. */
        for(iw=0;iw<nw;++iw)
        {
            x.push_back(d(band,iw,0));
            y.push_back(d(band,iw,1));
            z.push_back(d(band,iw,2));
        }
        /* These are required attributes from BasicTimeSeries.  We have
           to assume we can estract them from the components we just
           built. */
        ns=x[0].ns;
        tref=x[0].tref;
        t0=x[0].t0;
        live=false;   // Set this way in case we throw an error
        /* This is a sanity check.   Perhaps unnecessary overhead, but
         cost is not high for stability it can buy*/
        for(iw=0;iw<nw;++iw)
        {
            if((x[iw].ns!=ns) || (y[iw].ns!=ns) || (z[iw].ns!=ns) )
            {
                stringstream ss;
                ss << "Size mismatch in ComplexTimeSeries components for wavelet="
                    << iw<<endl
                    << "Testing against ns="<<ns<<endl
                    << "This set of wavelets has ns=("<<x[iw].ns
                    <<","<<y[iw].ns<<","<<z[iw].ns<<")"<<endl;;
                throw SeisppError(base_error + ss.str());
            }
            if((x[iw].t0!=t0) || (y[iw].t0!=t0) || (z[iw].t0!=t0) )
            {
                stringstream ss;
                ss << "Start time mismatch in ComplexTimeSeries components for wavelet="
                    << iw<<endl;
                throw SeisppError(base_error + ss.str());
            }
        }
        vector<ParticleMotionEllipse> pmi;   // nw estimates for each i
        ParticleMotionEllipse avg;
        ParticleMotionError err;
        pmi.reserve(nw);
        for(i=0;i<ns;++i)
        {
            for(iw=0;iw<nw;++iw)
            {
                pmi.push_back(ParticleMotionEllipse(x[iw].s[i],y[iw].s[i],z[iw].s[i],up));
            }
            ComputePMStats(pmi,avg,err);
            pmdata.push_back(avg);
            pmerr.push_back(err);
            pmi.clear();
        }
        this->post_attributes_to_metadata();
        live=true;
    }catch(...){throw;};
}
PMTimeSeries::PMTimeSeries(const PMTimeSeries& parent)
    : BasicTimeSeries(dynamic_cast<const BasicTimeSeries&> (parent)),
            Metadata(dynamic_cast<const Metadata&> (parent)),
              pmdata(parent.pmdata),  pmerr(parent.pmerr)
{
    f0=parent.f0;
    fw=parent.fw;
    averaging_length=parent.averaging_length;
    decfac=parent.decfac;
}

PMTimeSeries& PMTimeSeries::operator=(const PMTimeSeries& parent)
{
    if(this!=&parent)
    {
        f0=parent.f0;
        fw=parent.fw;
        averaging_length=parent.averaging_length;
        decfac=parent.decfac;
        this->BasicTimeSeries::operator=(parent);
        this->Metadata::operator=(parent);
        this->pmdata=parent.pmdata;
        this->pmerr=parent.pmerr;
    }
    return *this;
}


vector<ParticleMotionEllipse> PMTimeSeries::get_pmdata() {
    return(pmdata);
};

vector<ParticleMotionError> PMTimeSeries::get_pmerr() {
    return(pmerr);
};


/* Both of the following routines could use the at() method of std::vector
   but I use a custom test to allow me to throw a SeisppError, which is
   consistent with the rest of this code. */
ParticleMotionEllipse& PMTimeSeries::ellipse(int i)
{
    if(i<0 || i>=pmdata.size())
    {
        stringstream ss;
        ss << "PMTimeSeries::ellipse method:  "
            << "request for sample number "<<i
            << " is outside data range of "<< pmdata.size()<<endl;
        throw SeisppError(ss.str());
    }
    else
        return(pmdata[i]);
}
ParticleMotionError PMTimeSeries::errors(int i)
{
    if(i<0 || i>=pmerr.size())
    {
        stringstream ss;
        ss << "PMTimeSeries::errors method:  "
            << "request for sample number "<<i
            << " is outside data range of "<< pmerr.size()<<endl;
        throw SeisppError(ss.str());
    }
    else
        return(pmerr[i]);
}
/*! Helper procedure.  Wrapper function for C libmultiwavelet
  routine to estimate errors in particle motion ellipse parameters. 
  This procedure acts like a FORTRAN subroutine in that the average
  particle motion and error estimates are returned as arguments.  

arguments:
  d - ensemble of ParticleMotionEllipse objects to be averaged
  avg - average ParticleMotionEllipse
  err - errors
  */
void ComputePMStats(vector<ParticleMotionEllipse>& d, 
        ParticleMotionEllipse& avg, ParticleMotionError& err)
{
    /* This wrapper is a hideous inefficiency, but preferable to 
       rewriting pmvector_average which is quite complicated.  */
    int nd=d.size();
    int i,j;
    Particle_Motion_Ellipse *pmv;
    pmv = new Particle_Motion_Ellipse[nd];
    vector<double> major_amps, minor_amps;
    major_amps.reserve(nd);
    minor_amps.reserve(nd);
    for(i=0;i<nd;++i)
    {
        for(j=0;j<3;++j)
        {
            pmv[i].major[j]=d[i].major[j];
            pmv[i].minor[j]=d[i].minor[j];
        }
        pmv[i].rectilinearity=d[i].rectilinearity();
        major_amps.push_back(d[i].majornrm);
        minor_amps.push_back(d[i].minornrm);
    }
    Particle_Motion_Ellipse avgC;
    Particle_Motion_Error errC;
    pmvector_average(pmv,nd,&avgC,&errC);
    /* Now compute amplitude average and stdev as a 
       simple average of log values.  Error will be
       converted to dB at the end. */
    double sumlnamp,sumlndelta;
    for(i=0,sumlnamp=0.0;i<nd;++i)
        sumlnamp += log10(major_amps[i]);
    double major_amp_avg=sumlnamp/((double)nd);
    for(i=0,sumlndelta=0.0;i<nd;++i)
        sumlndelta = log10(major_amps[i]) / major_amp_avg;
    double major_err=sumlndelta/((double)(nd-1));
    /* Get rid of log for average */
    avg.majornrm=pow(10.0,major_amp_avg);
    /* Convert amplitude error to dB */
    err.dmajornrm=20.0*major_err;
    /* Now same for minor axis - repetitious but choice to not
       add function overhead*/
    for(i=0,sumlnamp=0.0;i<nd;++i)
        sumlnamp += log10(minor_amps[i]);
    double minor_amp_avg=sumlnamp/((double)nd);
    for(i=0,sumlndelta=0.0;i<nd;++i)
        sumlndelta = log10(minor_amps[i]) / minor_amp_avg;
    double minor_err=sumlndelta/((double)(nd=1));
    avg.minornrm=pow(10.0,minor_amp_avg);
    err.dminornrm=20.0*minor_err;
    /* Fill in particle motion average attributes */
    for(j=0;j<3;++j)
    {
        avg.major[j]=avgC.major[j];
        avg.minor[j]=avgC.minor[j];
    }
    /* Finally do the same for the error object. */
    err.dtheta_major=errC.dtheta_major;
    err.dphi_major=errC.dphi_major;
    err.dtheta_minor=errC.dtheta_minor;
    err.dphi_minor=errC.dphi_minor;
    err.delta_rect=errC.delta_rect;
    err.ndgf_major=errC.ndgf_major;
    err.ndgf_minor=errC.ndgf_minor;
    err.ndgf_rect=errC.ndgf_rect;
    err.ndgf_major_amp=nd-1;
    err.ndgf_minor_amp=nd-1;
}
void PMTimeSeries::zero_gaps()
{
	double tsend;
	int i,istart,iend;
	set<TimeWindow,TimeWindowCmp>::iterator this_gap;
	tsend = t0+((double)(ns-1))*dt;

	for(this_gap=gaps.begin();this_gap!=gaps.end();++this_gap)
	{
        	if(this_gap->end < t0) continue;
        	if(this_gap->start > tsend) continue;
        	if(this_gap->start<t0)
        		istart = 0;
        	else
        	    istart = SEISPP::nint((this_gap->start-t0)/dt);
        	if(this_gap->end>tsend)
        	    iend = ns-1;
        	else
     	            iend = SEISPP::nint((this_gap->end-t0)/dt);
        	for(i=0;i<3;++i)
        	    for(int j=istart;j<=iend;++j)
                {
                    pmdata[j].zero();
                    pmerr[j].zero();
                }
    	}
}
/* This is a series of methods to retrieve simplified representations of
 * the particle motion data as a function of time as time series 
 * objects.  In all cases the metadata is a copy of "this".   
 * This series of method are painfully paraellel.   There is probably
 * a way to have made this a template, but because the template 
 * variable is a method name I wasn't sure how to do that.*/
TimeSeries PMTimeSeries::rectilinearity()
{
    try{
        TimeSeries dts(dynamic_cast<Metadata&>(*this),false);
        dts.BasicTimeSeries::operator=(dynamic_cast<BasicTimeSeries&>(*this));
        /* set this so we can have a clue what this object is. 
         * PMDerivedTSType is defined in PMTimeSeries.h */
        dts.put(PMDerivedTSType,"rectilinearity");
        /* The following would normally be useful but the TimeSeries
         * constuctor does this so we drop it.  Left here in case 
         * that changes. 
        dts.reserve(this->ns);  
        */
        int i,k;
        for(i=0;i<(this->ns);++i)
        {
            ParticleMotionEllipse pmd=this->ellipse(i);
            dts.s.push_back(pmd.rectilinearity());
        }
        /* before returning be sure this is set true.   Original constructor
         * does not set this flag*/
        dts.live=true;
        return dts;
    }catch(...){throw;};
}
TimeSeries PMTimeSeries::major_axis_amplitude()
{
    try{
        TimeSeries dts(dynamic_cast<Metadata&>(*this),false);
        dts.BasicTimeSeries::operator=(dynamic_cast<BasicTimeSeries&>(*this));
        /* set this so we can have a clue what this object is. 
         * PMDerivedTSType is defined in PMTimeSeries.h */
        dts.put(PMDerivedTSType,"major_axis_amplitude");
        /* The following would normally be useful but the TimeSeries
         * constuctor does this so we drop it.  Left here in case 
         * that changes. 
        dts.reserve(this->ns);  
        */
        int i,k;
        for(i=0;i<(this->ns);++i)
        {
            ParticleMotionEllipse pmd=this->ellipse(i);
            dts.s.push_back(pmd.majornrm);
        }

        /* before returning be sure this is set true.   Original constructor
         * does not set this flag*/
        dts.live=true;
        return dts;
    }catch(...){throw;};
}
TimeSeries PMTimeSeries::minor_axis_amplitude()
{
    try{
        TimeSeries dts(dynamic_cast<Metadata&>(*this),false);
        dts.BasicTimeSeries::operator=(dynamic_cast<BasicTimeSeries&>(*this));
        /* set this so we can have a clue what this object is. 
         * PMDerivedTSType is defined in PMTimeSeries.h */
        dts.put(PMDerivedTSType,"minor_axis_amplitude");
        /* The following would normally be useful but the TimeSeries
         * constuctor does this so we drop it.  Left here in case 
         * that changes. 
        dts.reserve(this->ns);  
        */
        int i,k;
        for(i=0;i<(this->ns);++i)
        {
            ParticleMotionEllipse pmd=this->ellipse(i);
            dts.s.push_back(pmd.minornrm);
        }

        /* before returning be sure this is set true.   Original constructor
         * does not set this flag*/
        dts.live=true;
        return dts;
    }catch(...){throw;};
}
TimeSeries PMTimeSeries::major_azimuth()
{
    try{
        TimeSeries dts(dynamic_cast<Metadata&>(*this),false);
        dts.BasicTimeSeries::operator=(dynamic_cast<BasicTimeSeries&>(*this));
        /* set this so we can have a clue what this object is. 
         * PMDerivedTSType is defined in PMTimeSeries.h */
        dts.put(PMDerivedTSType,"major_azimuth");
        /* The following would normally be useful but the TimeSeries
         * constuctor does this so we drop it.  Left here in case 
         * that changes. 
        dts.reserve(this->ns);  
        */
        int i,k;
        for(i=0;i<(this->ns);++i)
        {
            ParticleMotionEllipse pmd=this->ellipse(i);
            dts.s.push_back(pmd.major_azimuth());
        }

        /* before returning be sure this is set true.   Original constructor
         * does not set this flag*/
        dts.live=true;
        return dts;
    }catch(...){throw;};
}
TimeSeries PMTimeSeries::major_inclination()
{
    try{
        TimeSeries dts(dynamic_cast<Metadata&>(*this),false);
        dts.BasicTimeSeries::operator=(dynamic_cast<BasicTimeSeries&>(*this));
        /* set this so we can have a clue what this object is. 
         * PMDerivedTSType is defined in PMTimeSeries.h */
        dts.put(PMDerivedTSType,"major_inclination");
        /* The following would normally be useful but the TimeSeries
         * constuctor does this so we drop it.  Left here in case 
         * that changes. 
        dts.reserve(this->ns);  
        */
        int i,k;
        for(i=0;i<(this->ns);++i)
        {
            ParticleMotionEllipse pmd=this->ellipse(i);
            dts.s.push_back(pmd.major_inclination());
        }

        /* before returning be sure this is set true.   Original constructor
         * does not set this flag*/
        dts.live=true;
        return dts;
    }catch(...){throw;};
}
TimeSeries PMTimeSeries::minor_azimuth()
{
    try{
        TimeSeries dts(dynamic_cast<Metadata&>(*this),false);
        dts.BasicTimeSeries::operator=(dynamic_cast<BasicTimeSeries&>(*this));
        /* set this so we can have a clue what this object is. 
         * PMDerivedTSType is defined in PMTimeSeries.h */
        dts.put(PMDerivedTSType,"minor_azimuth");
        /* The following would normally be useful but the TimeSeries
         * constuctor does this so we drop it.  Left here in case 
         * that changes. 
        dts.reserve(this->ns);  
        */
        int i,k;
        for(i=0;i<(this->ns);++i)
        {
            ParticleMotionEllipse pmd=this->ellipse(i);
            dts.s.push_back(pmd.minor_azimuth());
        }

        /* before returning be sure this is set true.   Original constructor
         * does not set this flag*/
        dts.live=true;
        return dts;
    }catch(...){throw;};
}
TimeSeries PMTimeSeries::minor_inclination()
{
    try{
        TimeSeries dts(dynamic_cast<Metadata&>(*this),false);
        dts.BasicTimeSeries::operator=(dynamic_cast<BasicTimeSeries&>(*this));
        /* set this so we can have a clue what this object is. 
         * PMDerivedTSType is defined in PMTimeSeries.h */
        dts.put(PMDerivedTSType,"minor_inclination");
        /* The following would normally be useful but the TimeSeries
         * constuctor does this so we drop it.  Left here in case 
         * that changes. 
        dts.reserve(this->ns);  
        */
        int i,k;
        for(i=0;i<(this->ns);++i)
        {
            ParticleMotionEllipse pmd=this->ellipse(i);
            dts.s.push_back(pmd.minor_inclination());
        }

        /* before returning be sure this is set true.   Original constructor
         * does not set this flag*/
        dts.live=true;
        return dts;
    }catch(...){throw;};
}
ostream& operator<<(ostream& os, PMTimeSeries& d)
{
    os << dynamic_cast<Metadata&>(d)<<endl;
    int i,k;
    for(i=0;i<d.ns;++i)
    {
        ParticleMotionEllipse pmd=d.ellipse(i);
        ParticleMotionError pme=d.errors(i);
        for(k=0;k<3;++k) os << pmd.major[k]*pmd.majornrm<<" ";
        for(k=0;k<3;++k) os << pmd.minor[k]*pmd.minornrm<<" ";
        os << pme.dtheta_major<<" "
            << pme.dphi_major <<" "
            << pme.dtheta_minor <<" "
            << pme.dphi_minor <<" "
            << pme.dmajornrm <<" "
            << pme.dminornrm <<" "
            << pme.delta_rect <<" "
            << pme.ndgf_major <<" "
            << pme.ndgf_minor <<" "
            << pme.ndgf_rect <<" "
            << pme.ndgf_major_amp <<" "
            << pme.ndgf_minor_amp <<endl;
    }
    return os;
}
