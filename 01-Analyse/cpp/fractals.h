#include "typedef.h"

#ifndef FRACTALS_H
#define FRACTALS_H

//Général
cplx frac_mandelbrot_set(cplx (*func)(cplx),cplx z,int kmax,int modmax);
cplx frac_julia_set     (cplx (*func)(cplx),cplx z,int kmax,int modmax,cplx c);

//f=z^2
cplx frac_mandelbrot    (cplx z,int kmax);
cplx frac_julia         (cplx z,int kmax,cplx c);


#endif // FRACTALS_H
