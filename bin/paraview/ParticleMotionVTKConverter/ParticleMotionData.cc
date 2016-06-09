#include <fstream>
#include "perf.h"
#include "ParticleMotionData.h"
ParticleMotionData::ParticleMotionData(ThreeComponentSeismogram& d,
        RegionalCoordinates& regcoord) 
       : ThreeComponentSeismogram(d), coords(regcoord)
{
    scale=1.0;  // default this to 1 - useless but necessary initialization
    try {
        double lat,lon,elev;
        /* fetch these assuming they are degrees then convert to radians*/
        lat=this->get_double("site.lat");
        lon=this->get_double("site.lon");
        elev=this->get_double("site.elev");
        x0_geo.lat=rad(lat);
        x0_geo.lon=rad(lon);
        x0_geo.r=r0_ellipse(x0_geo.lat)+elev;
        Cartesian_point cp=coords.cartesian(x0_geo);
        x0[0]=cp.x1;
        x0[1]=cp.x2;
        x0[2]=cp.x3;
        //DEBUG
        cout << "Station "<< this->get_string("sta")
            << " has geo coords: "<< lat<<", "<<lon<<", "<<elev<<endl
            << "Cartesian location="<<x0[0]<<", "<<x0[1]<<", "<<x0[2]<<endl;
        /*
        dmatrix utmp;
        utmp=tr(this->u);
        cout << "Data matrix in constructor"<<endl
            << utmp<<endl;
            */
    }catch(...){throw;};
}
ParticleMotionData::ParticleMotionData(const ParticleMotionData& parent)
    : ThreeComponentSeismogram(parent)
{
    scale=parent.scale;
    coords=parent.coords;
    x0_geo=parent.x0_geo;
    for(int i=0;i<3;++i) x0[i]=parent.x0[i];
}
ParticleMotionData& ParticleMotionData::operator=(const ParticleMotionData& parent)
{
    if(this!=&parent)
    {
        try{
            /* Trick to simplify copy of base data */
            this->ThreeComponentSeismogram::operator=(parent);
            scale=parent.scale;
            coords=parent.coords;
            x0_geo=parent.x0_geo;
            for(int i=0;i<3;++i) x0[i]=parent.x0[i];
        }catch(...){throw;};
    }
    return(*this);
}

/* This series of methods build from this base method */
dmatrix ParticleMotionData::raw_motion(TimeWindow tw)
{
    const string base_error("ParticleMotionData::raw_motion(TimeWindow):  ");
    int ss,se;
    ss=this->sample_number(tw.start);
    se=this->sample_number(tw.end);
    if(ss<0)
    {
        stringstream ss;
        ss << base_error << "Requested window start time="
            << tw.start<<" is before data start time of "
            << t0;
        throw SeisppError(ss.str());
    }
    else if(se>=u.columns())
    {
        stringstream ss;
        ss << base_error << "Requested window end = "
            << tw.end << " is beyond the end of data at time "
            << this->endtime();
        throw SeisppError(ss.str());
    }
    int ncol=se-ss+1;
    dmatrix result(3,ncol);
    int i,j;
    for(j=0;j<ncol;++j)
        for(i=0;i<3;++i)
            result(i,j)=this->u(i,ss+j);
    //DEBUG
    /*
    cout << "data in raw_motion method"<<endl;
    dmatrix temp=tr(result);
    cout << temp <<endl;
    */
    return(result);
}
dmatrix ParticleMotionData::raw_motion()
{
    TimeWindow tw;
    tw.start=this->t0;
    tw.end=this->endtime();
    return(this->raw_motion(tw));
}
/* helper procedure for converting raw to scaled data */
void ParticleMotionData::cook(dmatrix& d)
{
    int ncol=d.columns();
    int i,j;
    for(j=0;j<ncol;++j)
    {
        for(i=0;i<3;++i)
        {
            d(i,j)=scale*d(i,j)+x0[i];
        }
    }
}
dmatrix ParticleMotionData::particle_motion(TimeWindow tw)
{
    try{
        dmatrix result(this->raw_motion(tw));
        this->cook(result);
        return(result);
    }catch(...){throw;};
}
dmatrix ParticleMotionData::particle_motion()
{
    dmatrix result(this->raw_motion());
    this->cook(result);
    return(result);
}
double ParticleMotionData::max_amplitude()
{
    try {
        dmatrix rd(this->raw_motion());
        double dmax(0.0);
        int i;
        for(i=0;i<this->ns;++i)
        {
            double d=dnrm2(3,this->u.get_address(0,i),1);
            if(d>dmax) dmax=d;
        }
        return(dmax*scale);
    }catch(...){throw;};
}
double ParticleMotionData::max_amplitude(TimeWindow tw)
{
    try {
        dmatrix rd(this->raw_motion(tw));
        double dmax(0.0);
        int i;
        int nc=rd.columns();
        for(i=0;i<nc;++i)
        {
            double d=dnrm2(3,this->u.get_address(0,i),1);
            if(d>dmax) dmax=d;
        }
        return(dmax*scale);
    }catch(...){throw;};
}
