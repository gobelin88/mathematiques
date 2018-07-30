#include <math.h>
#include <complex>
#include <limits>
#include "Hypercplx.h"

#ifndef TYPEDEF_H
#define TYPEDEF_H

#define TYPE double
#define TYPE_EPSILON std::numeric_limits<TYPE>::epsilon()*10
#define TYPE_MAX std::numeric_limits<TYPE>::max()
#define TYPE_MIN std::numeric_limits<TYPE>::lowest()

#define LN2                 TYPE(0.6931471805599453094172321214581765680755)
#define EULAR_MACHERONIS    TYPE(0.5772156649015328606065120900824024310421)
#define M_PI                TYPE(3.1415926535897932384626433832795028841971)
#define M_PI2               TYPE(2*M_PI)
#define PI2_6               (M_PI*M_PI/6)

typedef std::complex<double> cplx;
//typedef HyperCplx<-1> cplx;

#define I cplx(0,1)

#endif // TYPEDEF_H

