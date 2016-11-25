#include <cfloat>
#include "ParticleMotionEllipse.h"
/* We need this prototypes from perf.h - BLAS functions for doubles.   Necessary
   to avoid nasty collisions with keyword complex when perf.h is incluced. */
extern "C" {
double ddot(int n,double *x, int incx, double *y, int incy);
void dscal(int n, double a, double *x, int incx);
double dnrm2(int n, double *x, int incx);
void dcopy(int n,double *x, int incx, double *y, int incy);
void daxpy(int n,double a,double *x, int incx, double *y, int incy);
void cgesvd ( char jobu, char jobvt, int m, int n, 
        FORTRAN_complex *ca, int lda, float *s, FORTRAN_complex *cu, int ldu, 
        FORTRAN_complex *cvt, int ldvt, int *info );
}
/* Default constructor forces initialization to 0.0.   */

ParticleMotionEllipse::ParticleMotionEllipse()
{
    int k;
    for(k=0;k<3;++k) major[k]=0.0;
    for(k=0;k<3;++k) minor[k]=0.0;
    majornrm=0.0;
    minornrm=0.0;
}
ParticleMotionEllipse::ParticleMotionEllipse(double *majin, double *minin)
{
    majornrm=dnrm2(3,majin,1);
    minornrm=dnrm2(3,minin,1);
    for(int i=0;i<3;++i)
    {
        major[i]=majin[i]/majornrm;
        minor[i]=minin[i]/minornrm;
    }
}

/*
   This constructor is a C++ conversion of a procedure in the older
libmultiwavelet library with the name compute_particle_motions.   
It uses an novel formula that works to compute a particle motion
ellipse from a narrow frequency signal using the phase between 
each of the 3 components x,y,z.   It thus manipulates 3 complex
numbers using formula found in this code.

This is the original documentation of the old procedure:

This function takes complex numbers x, y, and z defined by an 
eigenvector for a multiwavelet (would work for Fourier transforms
too, however) and returns a pointer to a structure that defines
the major and minor axes of the particle motion vectors defined
by those three complex numbers.  The up vector (assumed to be
three element vector) defines the direction used to resolve
the sign ambiguity inherent in defining an ellipse.  That is,
both the major and minor component directions are required
to have a positive projection in the up direction.  If they 
aren't the sign is flipped before returning.  Normally up 
would point [0,0,1] or in the up radial direction for P waves.
For S, it becomes more ambiguous and should be sorted out 
by a more complicated method.

The polarization information (defined by the ParticleMotionEllipse 
structure) is allocated within this routine. 

Author:  G. L. Pavlis
Written:  October 1999

Revision to C++ class constructor by pavlis, May 2015
*/
ParticleMotionEllipse::ParticleMotionEllipse(Complex x, Complex y,
        Complex z, double up[3])
{
	double rx,ry,rz,thetax,thetay,thetaz;  /* polar forms of x,y,z*/
	double a,b;
	double phi1,phi2;
	double x1[3],x2[3];
	double nrmx1,nrmx2,xmaj,xmin;


        rx=abs(x);
        ry=abs(y);
        rz=abs(z);
        /* Necessary for testing to avoid nans */
        if(rz<FLT_EPSILON && ry<FLT_EPSILON && rz<FLT_EPSILON )
        {
            this->zero();
            return;
        }

        thetax = atan2(std::imag(x),std::real(x));
        thetay = atan2(std::imag(y),std::real(y));
        thetaz = atan2(std::imag(z),std::real(z));

	a = rx*rx*cos(2.0*thetax) 
		+ ry*ry*cos(2.0*thetay) 
		+ rz*rz*cos(2.0*thetaz);
	b = rx*rx*sin(2.0*thetax) 
		+ ry*ry*sin(2.0*thetay) 
		+ rz*rz*sin(2.0*thetaz);

	phi1 = atan2(-b,a)/2.0;
	phi2 = phi1 + M_PI_2;

	x1[0] = rx*cos(phi1+thetax);
	x1[1] = ry*cos(phi1+thetay);
	x1[2] = rz*cos(phi1+thetaz);
	x2[0] = rx*cos(phi2+thetax);
	x2[1] = ry*cos(phi2+thetay);
	x2[2] = rz*cos(phi2+thetaz);

	nrmx1 = dnrm2(3,x1,1);
	nrmx2 = dnrm2(3,x2,1);
	/* normalize to unit vectors */
	dscal(3,1.0/nrmx1,x1,1);
	dscal(3,1.0/nrmx2,x2,1);

	if(nrmx1>nrmx2)
	{
		dcopy(3,x1,1,major,1);
		dcopy(3,x2,1,minor,1);
                majornrm=nrmx1;
                minornrm=nrmx2;
	}
	else
	{
		dcopy(3,x2,1,major,1);
		dcopy(3,x1,1,minor,1);
                majornrm=nrmx2;
                minornrm=nrmx1;
	}
	/* Choose the positive sign direction */
	if(ddot(3,up,1,major,1) < 0.0)
		dscal(3,-1.0,major,1);
	if(ddot(3,up,1,minor,1) < 0.0)
		dscal(3,-1.0,minor,1);
}
ParticleMotionEllipse::ParticleMotionEllipse(ComplexTimeSeries& x, 
        ComplexTimeSeries& y, 
            ComplexTimeSeries& z,
                TimeWindow w, double up[3])
{
    try {
        const string 
            base_error("ParticleMotionEllipse(ComplexTimeSeries contructor):  ");
        /* This is the work array used by cgesvd*/
        FORTRAN_complex *A,*U,*Vt;  // U and Vt are not used but required args 
        float svalues[3];
        /* Extract the windows*/
        ComplexTimeSeries xw=WindowData<ComplexTimeSeries>(x,w);
        ComplexTimeSeries yw=WindowData<ComplexTimeSeries>(y,w);
        ComplexTimeSeries zw=WindowData<ComplexTimeSeries>(z,w);
        int ntw=xw.s.size();  // assume all the same length
        A = (FORTRAN_complex *)calloc(3*ntw,sizeof(FORTRAN_complex));
        if(A==NULL) throw SeisppError(base_error
                + "calloc failed for principal component work matrix");
        /* Load A in fortran complex order and call the lapack svd routine*/
        int i,ia;
        for(i=0,ia=0;i<ntw;++i,ia+=3)
        {
            A[ia].r=(float)xw.s[i].real();
            A[ia].i=(float)xw.s[i].imag();
            A[ia+1].r=(float)yw.s[i].real();
            A[ia+1].i=(float)yw.s[i].imag();
            A[ia+2].r=(float)zw.s[i].real();
            A[ia+2].i=(float)zw.s[i].imag();
        }
        int info;
        cgesvd('o','n',3,ntw,A,3,svalues,U,3,Vt,3,&info);
        if(info!=0) 
        {
            free(A);
            throw SeisppError(base_error
                + "cgesvd returned an error");
        }
        /* We now use the first column of A that is overwritten
           with the singular vector linked to the largest singular
           value.  We scale the complex numbers by the singular
           value to get the amplitude */
        Complex xz,yz,zz;
        xz=Complex(A[0].r,A[0].i);
        yz=Complex(A[1].r,A[1].i);
        zz=Complex(A[2].r,A[2].i);
        /* This is very inefficient, but a simple way to build the 
           ellipse from the singular vector */
        ParticleMotionEllipse pmtmp(xz,yz,zz,up);
        for(i=0;i<3;++i)
        {
            this->major[i]=pmtmp.major[i];
            this->minor[i]=pmtmp.minor[i];
        }
        this->majornrm=pmtmp.majornrm;
        this->minornrm=pmtmp.minornrm;
        free(A);
    }catch(...){throw;};
}
/* Copy constructor */
ParticleMotionEllipse::ParticleMotionEllipse(const ParticleMotionEllipse& parent)
{
    int k;
    for(k=0;k<3;++k)
    {
        major[k]=parent.major[k];
        minor[k]=parent.minor[k];
    }
    majornrm=parent.majornrm;
    minornrm=parent.minornrm;
}
ParticleMotionEllipse& ParticleMotionEllipse::operator=(
        const ParticleMotionEllipse& parent)
{
    if(this!=&parent)
    {
        int k;
        for(k=0;k<3;++k)
        {
            major[k]=parent.major[k];
            minor[k]=parent.minor[k];
        }
        majornrm=parent.majornrm;
        minornrm=parent.minornrm;
    }
    return(*this);
}
double ParticleMotionEllipse::rectilinearity()
{
    if((minornrm>FLT_EPSILON) && (majornrm>FLT_EPSILON) )
        return(1.0 - minornrm/majornrm);
    else
        return(0.0);
}
double ParticleMotionEllipse::major_inclination()
{
    return(acos(major[2]));
}
double ParticleMotionEllipse::major_azimuth()
{
    return(M_PI_2 - atan2(major[1],major[0]));
}
double ParticleMotionEllipse::minor_inclination()
{
    return(acos(minor[2]));
}
double ParticleMotionEllipse::minor_azimuth()
{
    return(M_PI_2 - atan2(minor[1],minor[0]));
}
dmatrix ParticleMotionEllipse::points(int n)
{
    try {
        dmatrix result(n,3);
        /* x is a vector in principle coordinates and xp is
           used to hold the result after using a form of 
           transformation matrix implemented here with daxpy. */
        double xp[3];
        double dphi=2.0*M_PI/((double)n);
        int i,j;
        double phi;
        for(i=0,phi=0.0;i<n;++i,phi+=dphi)
        {
            for(j=0;j<3;++j) xp[j]=0.0;
            daxpy(3,majornrm*cos(phi),major,1,xp,1);
            daxpy(3,minornrm*sin(phi),minor,1,xp,1);
            for(j=0;j<3;++j) result(i,j)=xp[j];
        }
        return result;
    }catch(...){throw;};
}
void ParticleMotionEllipse::zero()
{
    for(int k=0;k<3;++k)
    {
        major[k]=0.0;
        minor[k]=0.0;
    }
    majornrm=0.0;
    minornrm=0.0;
}
ostream& operator<<(ostream& os, ParticleMotionEllipse& pmd)
{
    int k;
    for(k=0;k<3;++k)
        os << pmd.major[k]*pmd.majornrm<<" ";
    for(k=0;k<3;++k)
        os << pmd.minor[k]*pmd.minornrm<<" ";
    os <<endl;
    return os;
}
