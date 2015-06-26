#include "MWTransform.h"
#include "TimeSeries.h"
using namespace SEISPP;
MWTBundle::MWTBundle(ThreeComponentSeismogram& d,MWTransform& processor)
    : Metadata(dynamic_cast<Metadata&> (d))
{
    /* First make sure the transformation matrix is posted to metadata */
    put("U11",d.tmatrix[0][0]);
    put("U21",d.tmatrix[1][0]);
    put("U31",d.tmatrix[2][0]);
    put("U21",d.tmatrix[0][1]);
    put("U22",d.tmatrix[1][1]);
    put("U23",d.tmatrix[2][1]);
    put("U31",d.tmatrix[0][2]);
    put("U32",d.tmatrix[1][2]);
    put("U33",d.tmatrix[2][2]);
    /* The MW transform is a scalar operator so we need to split the 
       data into 3 scalar time series and then transform them.  */
    try {
        int i;
        this->mwtdata.reserve(3);
        for(i=0;i<3;++i)
        {
            TimeSeries *comp;
            try {
            comp=ExtractComponent(d,i);
            MWTMatrix mwt_this_comp(processor.transform(*comp));
            mwtdata.push_back(mwt_this_comp);
            } catch(...)
            {
                delete comp;
                throw;
            }
        }
        /* We can assume all the data in the matrix have the same
           number of bands and wavelets. */
        nw=processor.number_wavelet_pairs();
        nb=processor.number_frequencies();
    }catch(...){throw;};
}
MWTBundle::MWTBundle(TimeSeriesEnsemble& d,MWTransform& processor)
    : Metadata(dynamic_cast<Metadata&> (d))
{
    const string base_error("MWTBundle TimeSeriesEnsemble constructor:  ");
    /* We need some basic sanity checks on the ensemble */
    if(d.member.size()<=0) throw SeisppError(base_error
            + "No data - ensemble is empty.");
    int i;
    double dttest;
    bool has_live_members(false);
    vector<TimeSeries>::iterator dptr;
    for(dptr=d.member.begin(),i=0;dptr!=d.member.end();++dptr,++i)
    {
        if(dptr->live)
        {
            has_live_members=true;
            dttest=dptr->dt;
        }
        else if(has_live_members)
        {
            if(dptr->dt != dttest) 
                throw SeisppError(base_error
                        + "Sample rate mismatch in ensemble.\n"
                        + "All members must have a common sample rate");
        }
    }
    for(dptr=d.member.begin(),i=0;dptr!=d.member.end();++dptr,++i)
    {
        if(dptr->live)
        {
            try{
                mwtdata.push_back(processor.transform(*dptr));
            }catch(SeisppError& serr)
            {
                cerr << base_error << "MWTransform failed for member="<<i<<endl
                    << "Data from this member discarded."<<endl
                    << "Error message from MWTransform:"<<endl;
                serr.log_error();
            }
        }
    }
    if(mwtdata.size() <=0) throw SeisppError(base_error
            + "All ensemble members were either dead or failed processing");
    nw=processor.number_wavelet_pairs();
    nb=processor.number_frequencies();
}

MWTBundle::MWTBundle(const MWTBundle& parent)
    : mwtdata(parent.mwtdata)
{
    nw=parent.nw;
    nb=parent.nb;
}
MWTBundle& MWTBundle::operator=(const MWTBundle& parent)
{
    if(this!=&parent)
    {
        nw=parent.nw;
        nb=parent.nb;
        mwtdata=parent.mwtdata;
    }
    return(*this);
}

string MWTBundle::band_valid_test(int nbtest)
{
    if(nbtest>=0 && nbtest<nb)
        return string("ok");
    else
    {
        stringstream ss;
        ss << "Requested invalid band index="<<nbtest<<endl
            << "Must be between 0 and "<<nbtest-1<<endl;
        return(string(ss.str()));
    }
}
/* Now all the getters.  All assume member 0 is as good as any and they
also assume the ensemble is not empty. */
double MWTBundle::get_f0(int band)
{
    const string base_error("MWTBundle::get_f0:  ");
    string bok=band_valid_test(band);
    if(bok=="ok")
        return(mwtdata[0].get_f0(band));
    else
        throw SeisppError(base_error
                + bok);
}
double MWTBundle::get_fw(int band)
{
    const string base_error("MWTBundle::get_fw:  ");
    string bok=band_valid_test(band);
    if(bok=="ok")
        return(mwtdata[0].get_fw(band));
    else
        throw SeisppError(base_error
                + bok);
}
int MWTBundle::get_decfac(int band)
{
    const string base_error("MWTBundle::get_decfac:  ");
    string bok=band_valid_test(band);
    if(bok=="ok")
        return(mwtdata[0].get_decfac(band));
    else
        throw SeisppError(base_error
                + bok);
}
double MWTBundle::sample_interval(int band)
{
    const string base_error("MWTBundle::sample_interval:  ");
    string bok=band_valid_test(band);
    if(bok=="ok")
        return(mwtdata[0].sample_interval(band));
    else
        throw SeisppError(base_error
                + bok);
}
MWTwaveform MWTBundle::operator()(int b, int w, int m)
{
    const string base_error("MWTBundle::operator():");
    MWTMatrix mwtm;
    try{
        mwtm=mwtdata.at(m);
    }catch(std::exception& e)
    {
        stringstream ss;
        ss << "Request for data member="<<m<<" not consistent with "
            << "MWTBundle size="<<mwtdata.size();
        throw SeisppError(base_error+ss.str());
    }
    try {
        return mwtm(b,w);
    }catch(...){throw;};
}
Complex MWTBundle::operator()(int b, int w, int m, int iz)
{
    /* Not the most efficient way to implement this as this
       requires making a copy of a MWTwaveform, but this 
       routine is not a high priority at the time of initial
       development.  Could be changed if this becomes a 
       bottleneck.*/
    try {
        MWTwaveform d=this->operator()(b,w,m);
        return d.s[iz];
    }catch(...){throw;};
}

