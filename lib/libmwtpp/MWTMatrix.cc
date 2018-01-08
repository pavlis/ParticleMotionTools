#include "MWTransform.h"
MWTMatrix::MWTMatrix()
{
    nbands=0;
    nwavelets=0;
}
/* This is the main constructor for this object which is essentially
   an interface into the output of the existing multiwavelet transform
   C produre. */
MWTMatrix::MWTMatrix(MWtrace **draw,int nb, int nw,Metadata& md)
{
    this->Metadata::operator=(md);
    nbands=nb;
    nwavelets=nw;
    int ncomponents=nbands*nwavelets;
    d.reserve(ncomponents);
    /* Note **draw is nbands by nwavelets.  We preserve that concept
       but use FORTRAN order to define the matrix but using 
       implicit indexing through the private method matrix_index.
       Given how C does 2D arrays this effectively transposes the
       matrix in terms of how it is arranged in memory, but the
       interface makes this invisible.
     */
    int i,j;
    for(i=0;i<nbands;++i)
        for(j=0;j<nwavelets;++j)
        {
            MWTwaveform dij(draw[i][j]);
            d.push_back(dij);
        }
}
MWTMatrix::MWTMatrix(const MWTMatrix& parent) : Metadata(parent)
{
    nbands=parent.nbands;
    nwavelets=parent.nwavelets;
    d=parent.d;
}
vector<double> MWTMatrix::real(int band, int nw)
{
    string base_error("MWTMatrix::real:  ");
    string test=range_test(band,nw);
    if(test!="ok") 
        throw SeisppError(base_error+test);
    MWTwaveform work(d[matrix_index(band,nw)]);
    vector<double> result;
    result.reserve(work.ns);
    for(int i=0;i<work.ns;++i)
    {
        SEISPP::Complex val;
        val=work.s[i];
        result.push_back(val.real());
    }
    return result;
}
vector<double> MWTMatrix::imag(int band, int nw)
{
    string base_error("MWTMatrix::imag:  ");
    string test=range_test(band,nw);
    if(test!="ok") 
        throw SeisppError(base_error+test);
    MWTwaveform work(d[matrix_index(band,nw)]);
    vector<double> result;
    result.reserve(work.ns);
    for(int i=0;i<work.ns;++i)
    {
        SEISPP::Complex val;
        val=work.s[i];
        result.push_back(val.imag());
    }
    return result;
}
MWTwaveform& MWTMatrix::operator()(int nb,int nw)
{
    string base_error("MWTMatrix::operator():  ");
    string test=range_test(nb,nw);
    if(test!="ok") 
        throw SeisppError(base_error+test);
    return(d[nb*nwavelets+nw]);
}
MWTMatrix& MWTMatrix::operator=(const MWTMatrix& parent)
{
    if(this!=&parent)
    {
        this->Metadata::operator=(parent);
        nbands=parent.nbands;
        nwavelets=parent.nwavelets;
        d=parent.d;
    }
    return(*this);
}
double MWTMatrix::get_f0(int nb)
{
    if(nb>=0 && nb<nbands)
    {
        /* internal vector d is an implicit matrix stored with
         * multiwavelets in order first.  Hence, we need to compute
         * an offset with the nwavlets attribute to get the decimated
         * value of f0 for other than band 0 */
        int nbref=nb*nwavelets;
        return(d[nbref].get_f0());
    }
    else
    {
        stringstream ss;
        ss << "MWTMatrix::f0:  "
            << "Illegal request for band="<<nb
            << "Data range: number bands="<<nbands
            <<endl;
        throw SeisppError(ss.str());
    }
}
double MWTMatrix::get_fw(int nb)
{
    if(nb>=0 && nb<nbands)
    {
        /* as for gbet_f0 compute offset this way */
        int nbref=nb*nwavelets;
        return(d[nbref].get_fw());
    }
    else
    {
        stringstream ss;
        ss << "MWTMatrix::fw():  "
            << "Illegal request for band="<<nb
            << "Data range: number bands="<<nbands
            <<endl;
        throw SeisppError(ss.str());
    }
}
int MWTMatrix::get_decfac(int nb)
{
    if(nb>=0 && nb<nbands)
    {
        /* as for gbet_f0 compute offset this way */
        int nbref=nb*nwavelets;
        return(d[nbref].get_decfac());
    }
    else
    {
        stringstream ss;
        ss << "MWTMatrix::decfac():  "
            << "Illegal request for band="<<nb
            << "Data range: number bands="<<nbands
            <<endl;
        throw SeisppError(ss.str());
    }
}
double MWTMatrix::sample_interval(int nb)
{
    if(nb>=0 && nb<nbands)
    {
        int nbref=nb*nwavelets;
        double dt0=d[nbref].get_dt0();
        return(dt0*((double)d[nbref].get_decfac()));
    }
    else
    {
        stringstream ss;
        ss << "MWTMatrix::sample_interval(int band):  "
        << "Illegal request for band="<<nb
        << "Data range: number bands="<<nbands
        <<endl;
        throw SeisppError(ss.str());
    }
}
int MWTMatrix::get_wavelet_length(int nb)
{
    if(nb>=0 && nb<nbands)
    {
        return d[nb].get_wavelet_length();
    }
    else
    {
        stringstream ss;
        ss << "MWTMatrix::get_wavelet_length(int band):  "
        << "Illegal request for band="<<nb
        << "Data range: number bands="<<nbands
        <<endl;
        throw SeisppError(ss.str());
    }
}
string MWTMatrix::range_test(int ib, int iw)
{
    if( (ib<0) || (ib>=nbands) || (iw<0) || (iw>=nwavelets) )
    {
        stringstream ss;
        ss << "MWTMatrix index requested is out of range"<<endl
            << "Requested band="<<ib<<" and wavelet="<<iw<<endl
            << "Number of bands="<<nbands<<" and number of wavelets="
            << nwavelets;
        return(string(ss.str()));
    }
    else
        return(string("ok"));
}
