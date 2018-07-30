#include "typedef.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
#include <QDir>
#include <QImage>
#include <QColor>
#include <QPainter>
#include "Hypercplx.h"

typedef std::vector<TYPE (*)(TYPE)> RealFuncList;
TYPE clamp(TYPE value,TYPE min,TYPE max);
TYPE logabs(cplx s);

cplx to_riemann_sphere(cplx z);
cplx from_riemann_sphere(cplx z);
void sph_to_cart(cplx z,TYPE & X,TYPE & Y,TYPE & Z);

void fft(double *data, int nn, int isign);
void plot_spectrum(const char * file,
                   TYPE (*func)(TYPE),
                   TYPE res,
                   TYPE range_min,
                   TYPE range_max,
                   TYPE clamp_min=-1e100,
                   TYPE clamp_max=1e100);

void plot_real(const char * file,
               TYPE (*func)(TYPE),
               TYPE res,
               TYPE range_min,
               TYPE range_max,
               TYPE clamp_min=-1e100,
               TYPE clamp_max=1e100);

void plot_real(const char * file, const char *header,
               RealFuncList func,
               TYPE res,
               TYPE range_min,
               TYPE range_max,
               TYPE clamp_min=-1e100,
               TYPE clamp_max=1e100);

void plot(const char * file,
          cplx (*func)(cplx),
          TYPE (*op) (cplx),
          TYPE res,
          TYPE rangex_min,
          TYPE rangey_min,
          TYPE rangex_max,
          TYPE rangey_max,
          TYPE clamp_min=-1e100,
          TYPE clamp_max=1e100);

void plot(const char * file,
          cplx (*func)(cplx,int),
          TYPE (*op) (cplx),
          int it,
          TYPE res,
          TYPE rangex_min,
          TYPE rangey_min,
          TYPE rangex_max,
          TYPE rangey_max,
          TYPE clamp_min=-1e100,
          TYPE clamp_max=1e100);

void plot_it(const char * file,
             cplx (*func)(cplx,int),
             TYPE (*op) (cplx),
             int itmin,
             int itmax,
             int itstep,
             TYPE resx,TYPE resy,
             TYPE rangex_min,
             TYPE rangey_min,
             TYPE rangex_max,
             TYPE rangey_max,
             TYPE clamp_min=-1e100,
             TYPE clamp_max=1e100);

//IMAGES----------------------------------------------

enum ColorMode
{
    MODE_ARG_ONLY,
    MODE_MOD_ONLY,
    MODE_MOD_ISO_ONLY,
    MODE_ARG_AND_MOD_LIN,
    MODE_ARG_AND_MOD_ISO,
};

struct Resolution
{
    Resolution(unsigned int x){this->x=x;this->y=x;}
    Resolution(unsigned int x,unsigned int y){this->x=x;this->y=y;}
    unsigned int x,y;
};

void plot_riemann_csv(QString file,
                  std::function<cplx (cplx)> func,
                  Resolution res,
                  ColorMode mode=MODE_ARG_ONLY);

void plot_riemann(QString file,
                  std::function<cplx (cplx)> func,
                  Resolution res,
                  TYPE phi,TYPE theta,
                  ColorMode mode=MODE_ARG_ONLY);

void plot_img(QString file,
              std::function<cplx (cplx)> func,
              Resolution res,
              TYPE rangex_min,
              TYPE rangey_min,
              TYPE rangex_max,
              TYPE rangey_max,ColorMode mode=MODE_ARG_ONLY);

void plot_grid(QString file,
              std::function<cplx (cplx)> func,
              Resolution res,
               Resolution res_grid,
              TYPE rangex_min,
              TYPE rangey_min,
              TYPE rangex_max,
              TYPE rangey_max);

void plot_img_it(QString file,
                 cplx (*func)(cplx,int),
                 Resolution res,
                 int itmin,int n,
                 TYPE rangex_min,
                 TYPE rangey_min,
                 TYPE rangex_max,
                 TYPE rangey_max,ColorMode mode=MODE_ARG_ONLY);

QColor getColor(cplx value,TYPE & maxmod,TYPE & minmod,ColorMode mode);
void get_value_stat(cplx value_in,cplx & value_out,TYPE & maxmod,TYPE & minmod);
void build_image(QString file, Resolution res, cplx * values, TYPE & maxmod, TYPE & minmod, ColorMode mode);

void plot_compile(QString filename,QString dir,int nb_per_col);
