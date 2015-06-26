#include "MWTransform.h"
/* These come from the C file multiwavelet.h.   I do not want the whole
   multiwavlet.h file because it creates nasty symbol collisions.  Hence
   I only enter function prototypes referenced in this file.  This is a 
   DANGEROUS maintenance issue, so be aware.   note to self (glp) june 2015*/
extern "C"{
MWbasis *load_multiwavelets_pf(Pf *pf,int *nwavelets);
Tbl **build_decimation_objects(Tbl **filelists, int nbands, int *decfac);
Tbl **define_decimation(Pf *pf, int *nbands);
MWtrace **MWtransform(float *trace, double dt, double starttime, int nsamples,
               MWbasis *basis, int nbasis, Tbl **decimators, int nbands);
void free_MWtrace_matrix(MWtrace **t,int nrl, int nrh, int ncl,int nch);
}

MWTransform::MWTransform()
{
    nbands=0;
    nbasis=0;
    decimators=NULL;
    dec_fac=NULL;
}
MWTransform::MWTransform(string fname)
{
    const string base_error("MWTransform file constructor:  ");
    Pf *pf;
    pfread(const_cast<char *>(fname.c_str()),&pf);
    if(pf==NULL)
        throw SeisppError(base_error
                + "pfread failed for pf file="+fname);
    /* Note both basis_functions and nbasis are private object attributes.
       Note also that both of these two C functions will abort the program
       with elog_die if there are problems.  This is a bit harsh, but 
       the reality of this implementation.   May need to be tempered 
       for some applications */
    basis_functions=load_multiwavelets_pf(pf,&nbasis);
    /* The first Tbl ** is created by parsing the pf file.  It will
       be discarded after the build_decimation_objects C function is called
       immediately after */
    Tbl **decimator_definitions;
    decimator_definitions=define_decimation(pf,&nbands); //note this sets nbands
    dec_fac = new int[nbands];
    decimators=build_decimation_objects(decimator_definitions, nbands, dec_fac);
    /* We need to destroy the decimator_definitions Tbl ** because it is only
       an intermediatry. */
    int i;
    for(i=0;i<nbands;++i)
        freetbl(decimator_definitions[i],free);
    free(decimator_definitions);
}
MWTransform::MWTransform(const MWTransform& parent)
{
    const string error("MWTransform copy constructor:  is not implemented\n");
    throw SeisppError(error+"Coding error.   Call constructor instead of copying");
}
MWTransform::~MWTransform()
{
    delete [] dec_fac;
    int i;
    for(i=0;i<nbands;++i)
        freetbl(decimators[i],free);
    free(decimators);
    for(i=0;i<nbasis;++i)
    {
        free(basis_functions[i].r);
        free(basis_functions[i].i);
    }
    free(basis_functions);
}
MWTMatrix MWTransform::transform(TimeSeries& d)
{
    int i;
   /* The old C MWtransform routine uses a float array as data. We 
      have to first convert */
    float *fd;
    fd=new float[d.ns];
    for(i=0;i<d.ns;++i) fd[i]=(float)d[i];
    MWtrace **mwtraw;
    mwtraw=MWtransform(fd,d.dt,d.t0,d.ns,basis_functions,nbasis,decimators,nbands);
    delete [] fd;
    try {
        MWTMatrix result(mwtraw,nbands,nbasis,dynamic_cast<Metadata&>(d));
        free_MWtrace_matrix(mwtraw,0,nbands-1,0,nbasis-1);
        return result;
    }catch(...)
    {
        free_MWtrace_matrix(mwtraw,0,nbands-1,0,nbasis-1);
        throw;
    }
}
vector<Complex> MWTransform::basis(int n)
{
    if(n<0 || n>nbasis)
    {
        stringstream ss;
        ss << " MWTransform::basis method:  "
            << "Request for wavelet number "<<n
            <<" is illegal"<<endl
            << "Number of basis functions in this set="<<nbasis<<endl;
        throw SeisppError(ss.str());
    }
    vector<Complex> result;
    int nz;
    nz=basis_functions[n].n;
    result.reserve(nz);
    int i;
    for(i=0;i<nz;++i)
    {
        Complex val(basis_functions[n].r[i],basis_functions[n].i[i]);
        result.push_back(val);
    }
    return(result);
}
MWTransform& MWTransform::operator=(const MWTransform& parent)
{
    const string base_error("MWTransform assignment operator:  ");
    throw SeisppError(base_error + "not implemented\nCoding error");
}
