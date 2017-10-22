#include <math.h>
#include <cfloat>
#include "PMTimeSeries.h"
#include "Vector3DBootstrapError.h"
using namespace SEISPP;
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
/* Helper procedure - returns a vector of data that are input
   vector values converted to decibels. Throws an error if 
   a value is negative unless it is very tiny  - defines as 
   number less than FLT_EPSILON */
vector<double> dbamp(vector<double>& x)
{
  int nx=x.size();
  vector<double> result;
  result.reserve(nx);
  vector<double>::iterator xptr;
  for(xptr=x.begin();xptr!=x.end();++xptr)
  {
    double valdb;
    if((*xptr)<0.0)
    {
      if(fabs(*xptr)<FLT_EPSILON)
        valdb=20.0*log10(FLT_EPSILON);
      else
        throw SeisppError(string("dbamp procedure: negative values in put array are nonsense"));
    }
    valdb=20.0*log10(*xptr);
    result.push_back(valdb);
  }
  return result;
}


/*! Helper procedure.  Wrapper function for C libmultiwavelet
  routine to estimate errors in particle motion ellipse parameters.
  This procedure acts like a FORTRAN subroutine in that the average
  particle motion and error estimates are returned as arguments.

  This is a major revision from an earlier implemntation that
  used a routine called pmvector_average in the old C multiwavelet
  library.  Found that procedure produced errors that were too
  large due to incorrect handling of multiple estimates of angles.
  This version uses a new bootstrap error approach.

arguments:
  d - ensemble of ParticleMotionEllipse objects to be averaged
  avg - average ParticleMotionEllipse
  err - errors
  */
void ComputePMStats(vector<ParticleMotionEllipse>& d,
        ParticleMotionEllipse& avg, ParticleMotionError& err,
        double confidence_level, int number_of_trials)
{
  try{
    /* This wrapper is a hideous inefficiency, but preferable to
       rewriting pmvector_average which is quite complicated.  */
    int nd=d.size();
    int i,j;
    vector<double> major_amps, minor_amps,rect;
    major_amps.reserve(nd);
    minor_amps.reserve(nd);
    rect.reserve(nd);
    /* These two hold normalized major and minor axis vectors */
    dmatrix dmajor(3,nd);
    dmatrix dminor(3,nd);
    for(i=0;i<nd;++i)
    {
      /* We assume the vectors passed are already normalized to be
      unit vectors - amplitude is contained in the majornrm and minornrm
      vectors */
        for(j=0;j<3;++j)
        {
            dmajor(j,i)=d[i].major[j];
            dminor(j,i)=d[i].minor[j];
        }
        major_amps.push_back(d[i].majornrm);
        minor_amps.push_back(d[i].minornrm);
        rect.push_back(d[i].rectilinearity());
    }
    /* This is the old procedure that computed errors.  Replacing here
    by bootstrap error estimation
    Particle_Motion_Ellipse avgC;
    Particle_Motion_Error errC;
    pmvector_average(pmv,nd,&avgC,&errC); */
    pair<double,double> majorampstats, minorampstats,rectstats;
    vector<double> xdb=dbamp(major_amps);
    majorampstats=bootstrap_mv(xdb,confidence_level,number_of_trials);
    avg.majornrm=majorampstats.first;
    err.dmajornrm=majorampstats.second;
    xdb=dbamp(minor_amps);
    minorampstats=bootstrap_mv(xdb,confidence_level,number_of_trials);
    avg.minornrm=minorampstats.first;
    err.dminornrm=minorampstats.second;
    rectstats=bootstrap_mv(rect,confidence_level,number_of_trials);
    //avg.rectilinearity=rectstats.first;
    err.delta_rect=rectstats.second;
    Vector3DBootstrapError majboot(dmajor,confidence_level,number_of_trials);
    Vector3DBootstrapError minboot(dminor,confidence_level,number_of_trials);
    /* For major we just copy the bootstrap mean */
    vector<double> vmed;
    vmed=majboot.mean_vector();
    for(j=0;j<3;++j) avg.major[j] = vmed[j];  // assumes vmed is unit vector
    /* For minor we want to force the vector to be perpendicular to major
    axis vector.  Borrowed form old code */
    double w[3],w2[3];
    dr3cros(avg.major,&(vmed[0]),w);
    dr3cros(w,avg.major,w2);
    for(j=0;j<3;++j) avg.minor[j]=w2[j];
    /* For this implementation we use the bootstrap error in the dot
    product angle between resampled observations as estimate for the
    error all angle terms.   Note these are retained as radians. */
    double aerr=majboot.angle_error();
    err.dtheta_major=aerr;
    err.dphi_major=aerr;
    aerr=majboot.angle_error();
    err.dtheta_minor=aerr;
    err.dphi_minor=aerr;
    /* For degrees of freedom we set all to nd-1 */
    err.ndgf_major=nd-1;
    err.ndgf_minor=nd-1;
    err.ndgf_rect=nd-1;
    err.ndgf_major_amp=nd-1;
    err.ndgf_minor_amp=nd-1;
  }catch(...){throw;};
}
PMTimeSeries::PMTimeSeries() : Metadata(), BasicTimeSeries()
{
    averaging_length=0;
    f0=0.0;
    fw=0.0;
    decfac=0;
}

PMTimeSeries::PMTimeSeries(MWTBundle& d, int band, int timesteps, int avlen,
    double confidence, int bsmultiplier)
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
        int nw=d.number_wavelets();
        int ntrials=bsmultiplier*nw;
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
                ComputePMStats(pmi,avg,err,confidence,ntrials);
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
PMTimeSeries::PMTimeSeries(MWTBundle& d, int band, double confidence,
        int bsmultiplier)
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
        int nw=d.number_wavelets();
        int ntrials=nw*bsmultiplier;
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
            ComputePMStats(pmi,avg,err,confidence,ntrials);
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
