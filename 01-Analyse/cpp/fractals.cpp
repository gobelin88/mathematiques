#include "fractals.h"

cplx frac_mandelbrot_set(cplx (*func)(cplx), cplx z, int kmax,int modmax)
{
    cplx s=z;
    int k=0;
    while(abs(s)<modmax && k++<kmax)
    {
        s=func(s)+z;
    }

    if(k>=kmax)return cplx(0,0);
    double theta=k/20.0;
    return cplx(sin(theta)*k,cos(theta)*k);
}

cplx frac_julia_set(cplx (*func)(cplx), cplx z, int kmax,int modmax,cplx c)
{
    cplx s=z;
    int k=0;
    while(abs(s)<modmax && k++<kmax)
    {
        s=func(s)+c;
    }

    if(k>=kmax)return cplx(0,0);
    double theta=k/20.0;
    return cplx(sin(theta)*k,cos(theta)*k);
}

cplx frac_mandelbrot(cplx z, int kmax)
{
    TYPE mz=sqrt((z.real()-0.25)*(z.real()-0.25)+z.imag()*z.imag());
    if(z.real()<(mz-2*mz*mz+0.25))return cplx(0,0);
    if(((z.real()+1)*(z.real()+1)+z.imag()*z.imag())<(1.0/16.0))return cplx(0,0);

    return frac_mandelbrot_set([](cplx s)->cplx{return s*s;},z,kmax,2);
}

cplx frac_julia(cplx z, int kmax,cplx c)
{
    return frac_julia_set([](cplx s)->cplx{return s*s;},z,kmax,2,c);
}
