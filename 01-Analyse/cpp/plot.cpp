#include "plot.h"

cplx to_riemann_sphere(cplx z)
{
    TYPE r=abs(z);
    TYPE theta=arg(z);
    return cplx( 2.0*atan(1.0/r) , theta );
}

cplx from_riemann_sphere(cplx z)
{
    TYPE r= sin(z.real())/(1.0-cos(z.real()));
    TYPE theta=z.imag();
    return cplx( r*sin(theta),r*cos(theta) );
}

void sph_to_cart (cplx z,TYPE & X,TYPE & Y,TYPE & Z)
{
    X=cos(z.imag())*sin(z.real());
    Y=sin(z.imag())*sin(z.real());
    Z=cos(z.real());
}

TYPE dot(TYPE X,TYPE Y,TYPE Z,TYPE vx,TYPE vy,TYPE vz)
{
    return vx*X+vy*Y+vz*Z;
}

void crossn(TYPE  ux,TYPE  uy,TYPE  uz,
            TYPE  vx,TYPE  vy,TYPE  vz,
            TYPE& wx,TYPE& wy,TYPE& wz)
{
    wx=uy*vz-uz*vy;
    wy=uz*vx-ux*vz;
    wz=ux*vy-uy*vx;

    TYPE n=sqrt(wx*wx+wy*wy+wz*wz);
    wx/=n;
    wy/=n;
    wz/=n;
}

TYPE clamp(TYPE value,TYPE min,TYPE max)
{
    if(std::isnan(value))return 0;
    else
    {
        if(value>max)return max;
        if(value<min)return min;
    }
    return value;
}

TYPE logabs(cplx s)
{
    return log(abs(s));
}
////////////////////////////////////////////////////////////

void fft(double * data, int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            tempr = data[j];     data[j] = data[i];     data[i] = tempr;
            tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n > mmax) {
        istep = 2*mmax;
        theta = M_PI2/(isign*mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j =i + mmax;
                tempr = wr*data[j]   - wi*data[j+1];
                tempi = wr*data[j+1] + wi*data[j];
                data[j]   = data[i]   - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr)*wpr - wi*wpi + wr;
            wi = wi*wpr + wtemp*wpi + wi;
        }
        mmax = istep;
    }
}

void plot_spectrum(const char * file,
                   TYPE (*func)(TYPE),
                   TYPE res,
                   TYPE range_min,
                   TYPE range_max,
                   TYPE clamp_min,
                   TYPE clamp_max)
{
    std::cout<<file<<std::endl;

    TYPE r=(int)pow(2.0, ceil(log((double)res)/log(2.0)));
    std::vector<double> sig(2*r),spec(2,r);
    for(int i=0;i<r;i++)
    {
        double x((range_max-range_min)*i/r+range_min);

        sig[2*i]=(clamp(func(x),clamp_min,clamp_max));
        sig[2*i+1]=0;

        std::cout<<i/(r)*100<<"\r";
        std::cout.flush();
    }
    spec=sig;
    fft(spec.data(),r,1);

    std::ofstream os(file);
    os<<"x;f(x);|F(x)|;arg(F(x))"<<std::endl;
    for(int i=0;i<r;i++)
    {
        double x((range_max-range_min)*i/r+range_min);
        os<<x<<";";
        os<<sig[2*i]<<";";
        os<<20*log10(abs(cplx(spec[2*i],spec[2*i+1])))<<";";
        os<<arg(cplx(spec[2*i],spec[2*i+1]))<<";"<<std::endl;
    }
    os.close();
}

void plot_real(const char * file,
               TYPE (*func)(TYPE),
               TYPE res,
               TYPE range_min,
               TYPE range_max,
               TYPE clamp_min,
               TYPE clamp_max)
{
    std::cout<<file<<std::endl;
    std::ofstream os(file);
    os<<"x;f(x)"<<std::endl;
    TYPE r=res;
    for(int i=0;i<=r;i++)
    {
        double x((range_max-range_min)*i/r+range_min);

        os<<x<<";"
         <<clamp(func(x),clamp_min,clamp_max)<<";"
        <<std::endl;

        std::cout<<i/(res)*100<<"\r";
        std::cout.flush();
    }
    os.close();
}

void plot_real(const char * file,
               const char * header,
               RealFuncList func,
               TYPE res,
               TYPE range_min,
               TYPE range_max,
               TYPE clamp_min,
               TYPE clamp_max)
{
    std::cout<<file<<std::endl;
    std::ofstream os(file);
    os<<header<<"\n";

    TYPE r=res;
    for(int i=0;i<=r;i++)
    {
        double x((range_max-range_min)*i/r+range_min);

        os<<x<<";";
        for(int j=0;j<func.size();j++)
        {
            os<<clamp(func[j](x),clamp_min,clamp_max)<<";";
        }
        os<<std::endl;

        std::cout<<i/(res)*100<<"\r";
        std::cout.flush();
    }
    os.close();
}

void plot(const char * file,
          cplx (*func)(cplx),
          TYPE (*op) (cplx),
          TYPE res,
          TYPE rangex_min,
          TYPE rangey_min,
          TYPE rangex_max,
          TYPE rangey_max,
          TYPE clamp_min,
          TYPE clamp_max)
{
    std::cout<<file<<std::endl;
    std::ofstream os(file);
    os<<"x;y;f(z)"<<std::endl;
    TYPE r=res;
    for(int i=0;i<=r;i++)
    {
        for(int j=0;j<=r;j++)
        {
            cplx z((rangex_max-rangex_min)*i/r+rangex_min,
                   (rangey_max-rangey_min)*j/r+rangey_min);

            os<<z.real()<<";"<<z.imag()<<";"
             <<clamp(op(func(z)),clamp_min,clamp_max)<<";"
            <<std::endl;


        }
        std::cout<<i/(res)*100<<"\r";
        std::cout.flush();
    }
    os.close();
}


void plot(const char * file,
          cplx (*func)(cplx,int),
          TYPE (*op) (cplx),
          int it,
          TYPE res,
          TYPE rangex_min,
          TYPE rangey_min,
          TYPE rangex_max,
          TYPE rangey_max,
          TYPE clamp_min,
          TYPE clamp_max)
{
    std::cout<<file<<std::endl;
    std::ofstream os(file);
    os<<"x;y;f(z)"<<std::endl;
    TYPE r=res;
    for(int i=0;i<=r;i++)
    {
        for(int j=0;j<=r;j++)
        {
            cplx z((rangex_max-rangex_min)*i/r+rangex_min,
                   (rangey_max-rangey_min)*j/r+rangey_min);

            os<<z.real()<<";"<<z.imag()<<";"
             <<clamp(op(func(z,it)),clamp_min,clamp_max)<<";"
            <<std::endl;


        }
        std::cout<<i/(res)*100<<"\r";
        std::cout.flush();
    }
    os.close();
}


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
             TYPE clamp_min,
             TYPE clamp_max)
{
    std::ofstream os(file);
    os<<"x;y;";
    for(int it=itmin;it<=itmax;it*=itstep)
    {
        os<<"it"<<it<<";";
    }
    os<<std::endl;

    double stepx=(rangex_max-rangex_min)/(resx);
    double stepy=(rangey_max-rangey_min)/(resy);

    std::cout<<file<<std::endl;

    for(int i=0;i<=resx;i++)
    {
        for(int j=0;j<=resy;j++)
        {
            cplx z(stepx*i+rangex_min,
                   stepy*j+rangey_min);

            os<<z.real()<<";"<<z.imag()<<";";

            for(int it=itmin;it<=itmax;it*=itstep)
            {
                os<<clamp(op(func(z,it)),clamp_min,clamp_max)<<";";
            }
            os<<std::endl;
        }
        std::cout<<i/(resx)*100<<"\r";
        std::cout.flush();
    }
    os.close();
}

//------------------------------------------------------------------------------
//
void plot_grid(QString file,
               std::function<cplx (cplx)> func,
               Resolution res,
               Resolution res_grid,
               TYPE rangex_min,
               TYPE rangey_min,
               TYPE rangex_max,
               TYPE rangey_max)
{
    QImage image(res.x,res.y,QImage::Format_RGB32);
    image.fill(QColor(255,255,255));
    QPainter painter(&image);

    std::cout<<file.toStdString().data()<<std::endl;

    painter.setPen(QColor(0,0,0));

    double dx=rangex_max-rangex_min;
    double dy=rangey_max-rangey_min;

    for(unsigned int i=0;i<res_grid.x;i++)
    {
        for(unsigned int j=0;j<res_grid.y;j++)
        {
            cplx z((dx*i)/res_grid.x+rangex_min,
                   (dy*j)/res_grid.y+rangey_min);
            cplx value=func(z);

            cplx z_img( ((z.real()-rangex_min)/dx)*res.x,
                        ((z.imag()-rangey_min)/dy)*res.y);
            cplx value_img( ((value.real()-rangex_min)/dx)*res.x,
                            ((value.imag()-rangey_min)/dy)*res.y);

            QPointF p1(z_img.real(),z_img.imag());
            QPointF p2(value_img.real(),value_img.imag());
            QPointF p3=p2-p1;

            double n3=sqrt(p3.x()*p3.x()+p3.y()*p3.y());

            if(n3>0)
            {
                double stepx=((double)res.x)/((double)res_grid.x);
                double stepy=((double)res.y)/((double)res_grid.y);
                double step=sqrt(stepx*stepx+stepy*stepy);

                painter.drawLine(p1-p3/(3*n3)*step,
                                 p1+p3/(3*n3)*step);

                painter.drawEllipse(p1-p3/(3*n3)*step,3,3);
            }
            else
            {
                painter.drawEllipse(p1,3,3);
            }

        }
        std::cout<<(i*100.0)/(res_grid.x)<<"\r";
        std::cout.flush();
    }

    image.save(file);
}

//------------------------------------------------------------------------------
//Img
void plot_riemann_csv(QString file,
                      std::function<cplx (cplx)> func,
                      Resolution res,
                      ColorMode mode)
{
    TYPE rangex_min=0.0;
    TYPE rangey_min=0.0;
    TYPE rangex_max=M_PI;
    TYPE rangey_max=2*M_PI;

    TYPE maxmod=std::numeric_limits<TYPE>::min(),minmod=std::numeric_limits<TYPE>::max();
    cplx * values=new cplx[res.x*res.y];

    std::cout<<file.toStdString().data()<<std::endl;
    std::ofstream os(file.toStdString());
    os<<"x;y;X;Y;Z;|f(z)|;f(z).arg;R;V;B"<<std::endl;
    for(unsigned int i=0;i<res.x;i++)
    {
        for(unsigned int j=0;j<res.y;j++)
        {
            int id=i+res.x*j;
            cplx zr((rangex_max-rangex_min)*i/res.x+rangex_min,
                    (rangey_max-rangey_min)*j/res.y+rangey_min);

            get_value_stat(func(from_riemann_sphere(zr)),values[id],maxmod,minmod);
        }
        std::cout<<i*100.0/(res.x)<<"\r";
        std::cout.flush();
    }

    for(unsigned int i=0;i<res.x;i++)
    {
        for(unsigned int j=0;j<res.y;j++)
        {
            int id=i+res.x*j;
            cplx zr((rangex_max-rangex_min)*i/res.x+rangex_min,
                    (rangey_max-rangey_min)*j/res.y+rangey_min);

            cplx z=from_riemann_sphere(zr);
            TYPE X,Y,Z;
            sph_to_cart(zr,X,Y,Z);

            QColor color=getColor(values[id],maxmod,minmod,mode);

            os<<z.real()<<";"<<z.imag()<<";";
            os<<X<<";"<<Y<<";"<<Z<<";";
            os<<values[id].real()<<";"<<values[id].imag()<<";";
            os<<color.red()<<";"<<color.green()<<";"<<color.blue()<<";";
            os<<std::endl;
        }
        std::cout<<i*100.0/(res.x)<<"\r";
        std::cout.flush();
    }

    os.close();

    delete [] values;
}

void plot_riemann(QString file,
                  std::function<cplx (cplx)> func,
                  Resolution res,
                  TYPE phi,TYPE theta,
                  ColorMode mode)
{
    TYPE rangex_min=0.0;
    TYPE rangey_min=0.0;
    TYPE rangex_max=M_PI;
    TYPE rangey_max=2*M_PI;

    TYPE maxmod=std::numeric_limits<TYPE>::min(),minmod=std::numeric_limits<TYPE>::max();
    cplx * values=new cplx[res.x*res.y];

    std::cout<<file.toStdString().data()<<std::endl;

    for(unsigned int i=0;i<res.x;i++)
    {
        for(unsigned int j=0;j<res.y;j++)
        {
            int id=i+res.x*j;
            cplx zr((rangex_max-rangex_min)*i/res.x+rangex_min,
                    (rangey_max-rangey_min)*j/res.y+rangey_min);

            get_value_stat(func(from_riemann_sphere(zr)),values[id],maxmod,minmod);
        }
        std::cout<<i*100.0/(res.x)<<"\r";
        std::cout.flush();
    }

    TYPE nx,ny,nz;
    TYPE vx,vy,vz;
    TYPE ux,uy,uz;
    sph_to_cart(cplx(phi,theta),nx,ny,nz);

    crossn(0,0,-1,
           nx,ny,nz,
           ux,uy,uz);

    crossn(ux,uy,uz,
           nx,ny,nz,
           vx,vy,vz);

    QImage image(res.x,res.y,QImage::Format_RGB32);
    image.fill(QColor(255,255,255));

    TYPE m=50;

    TYPE AXx=dot(2,0,0,ux,uy,uz);
    TYPE AXy=dot(2,0,0,vx,vy,vz);
    TYPE AYx=dot(0,2,0,ux,uy,uz);
    TYPE AYy=dot(0,2,0,vx,vy,vz);
    TYPE AZx=dot(0,0,2,ux,uy,uz);
    TYPE AZy=dot(0,0,2,vx,vy,vz);

    TYPE AXx1=dot(1,0,0,ux,uy,uz);
    TYPE AXy1=dot(1,0,0,vx,vy,vz);
    TYPE AYx1=dot(0,1,0,ux,uy,uz);
    TYPE AYy1=dot(0,1,0,vx,vy,vz);
    TYPE AZx1=dot(0,0,1,ux,uy,uz);
    TYPE AZy1=dot(0,0,1,vx,vy,vz);

    QPainter painter(&image);
    painter.setPen(QColor(0,0,0));
    painter.drawLine(QPointF((res.x-m*2)*0.5+m,res.y-2*m-1-(res.y-m*2)*0.5+m), QPointF((AXx+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(AXy+1)*(res.y-m*2)*0.5+m));
    painter.drawLine(QPointF((res.x-m*2)*0.5+m,res.y-2*m-1-(res.y-m*2)*0.5+m), QPointF((AYx+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(AYy+1)*(res.y-m*2)*0.5+m));
    painter.drawLine(QPointF((res.x-m*2)*0.5+m,res.y-2*m-1-(res.y-m*2)*0.5+m), QPointF((AZx+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(AZy+1)*(res.y-m*2)*0.5+m));

    for(unsigned int i=0;i<res.x;i++)
    {
        for(unsigned int j=0;j<res.y;j++)
        {
            int id=i+res.x*j;
            cplx zr((rangex_max-rangex_min)*i/res.x+rangex_min,
                    (rangey_max-rangey_min)*j/res.y+rangey_min);

            TYPE X,Y,Z;
            sph_to_cart(zr,X,Y,Z);

            TYPE dz=dot(X,Y,Z,nx,ny,nz);
            TYPE dx=dot(X,Y,Z,ux,uy,uz);
            TYPE dy=dot(X,Y,Z,vx,vy,vz);

            if(dz>0)
            {
                QPoint pix((dx+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(dy+1)*(res.y-m*2)*0.5+m);
                QRgb color=getColor(values[id],maxmod,minmod,mode).rgb();
                image.setPixel(pix,color);
                image.setPixel(pix.x()+1,pix.y(),color);
                image.setPixel(pix.x()-1,pix.y(),color);
                image.setPixel(pix.x(),pix.y()+1,color);
                image.setPixel(pix.x(),pix.y()-1,color);
                image.setPixel(pix.x()+1,pix.y()+1,color);
                image.setPixel(pix.x()-1,pix.y()-1,color);
                image.setPixel(pix.x()-1,pix.y()+1,color);
                image.setPixel(pix.x()+1,pix.y()-1,color);
            }
        }
        std::cout<<i*100.0/(res.x)<<"\r";
        std::cout.flush();
    }

    painter.drawLine(QPointF((AXx1+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(AXy1+1)*(res.y-m*2)*0.5+m), QPointF((AXx+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(AXy+1)*(res.y-m*2)*0.5+m));
    painter.drawLine(QPointF((AYx1+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(AYy1+1)*(res.y-m*2)*0.5+m), QPointF((AYx+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(AYy+1)*(res.y-m*2)*0.5+m));
    painter.drawLine(QPointF((AZx1+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(AZy1+1)*(res.y-m*2)*0.5+m), QPointF((AZx+1)*(res.x-m*2)*0.5+m,res.y-2*m-1-(AZy+1)*(res.y-m*2)*0.5+m));

    image.save(file);

    delete [] values;
}

void plot_img(QString file,
              std::function<cplx (cplx)> func,
              Resolution res,
              TYPE rangex_min,
              TYPE rangey_min,
              TYPE rangex_max,
              TYPE rangey_max,
              ColorMode mode)
{
    TYPE maxmod=std::numeric_limits<TYPE>::min(),minmod=std::numeric_limits<TYPE>::max();

    cplx * values=new cplx[res.x*res.y];

    std::cout<<file.toStdString().data()<<std::endl;
    for(unsigned int i=0;i<res.x;i++)
    {
        for(unsigned int j=0;j<res.y;j++)
        {
            int id=i+res.x*j;
            cplx z((rangex_max-rangex_min)*i/res.x+rangex_min,
                   (rangey_max-rangey_min)*j/res.y+rangey_min);

            cplx value=func(z);
            get_value_stat(value,values[id],maxmod,minmod);
        }
        std::cout<<(i*100.0)/(res.x)<<"\r";
        std::cout.flush();
    }

    build_image(file,res,values,maxmod,minmod,mode);

    delete [] values;
}

void plot_img_it(QString file,
                 cplx (*func)(cplx,int),
                 Resolution res,
                 int itmin,int n,
                 TYPE rangex_min,
                 TYPE rangey_min,
                 TYPE rangex_max,
                 TYPE rangey_max,
                 ColorMode mode)
{

    for(int it=itmin;it<=n;it*=2)
    {
        plot_img(QString(file).arg(it),[=](cplx z){return func(z,it);},res,rangex_min,rangey_min,rangex_max,rangey_max,mode);
    }
}

void get_value_stat(cplx value_in, cplx & value_out, TYPE & maxmod, TYPE & minmod)
{
    if(std::isnormal(value_in.real()))
    {
        value_out=cplx(log(1.0+abs(value_in)),arg(value_in));
        if(value_out.real()>maxmod)maxmod=value_out.real();
        if(value_out.real()<minmod)minmod=value_out.real();
    }
    else if(!std::isnan(value_in.real()))
    {
        value_out=cplx(log(1.0+abs(value_in)),arg(value_in));
    }
    else
    {
        value_out=NAN;
    }
}

//int extrapol(TYPE mod,TYPE minmod,TYPE moymod,TYPE maxmod,ColorMode mode)
//{
//    if(maxmod==minmod)return mod/minmod*128;
//    if(std::isnan(mod))return 0;

//    TYPE hmin=0;
//    TYPE hmoy=128;
//    TYPE hmax=255;
//    if(mode==MODE_FUNCTION)
//    {
//        hmin=64;
//    }
//    else if(mode==MODE_FRACTAL)
//    {
//        hmin=0;
//        moymod=(minmod+maxmod)/2.0;
//    }

//    TYPE m=0;
//    TYPE a=minmod;
//    TYPE b=moymod-m*(moymod-minmod);
//    TYPE c=moymod+m*(maxmod-moymod);
//    TYPE d=maxmod;

//    if(mod>=a && mod<b)
//    {
//        return (mod-a)/(b-a)*(hmoy-hmin)+hmin;
//    }
//    else if(mod>=b && mod<c)
//    {
//        return hmoy;
//    }
//    else if(mod>=c && mod<=d)
//    {
//        return (mod-c)/(d-c)*(hmax-hmoy)+hmoy;
//    }
//    else
//    {
//        return 0;
//    }
//}

//TYPE get_median(complexd * values,TYPE res)
//{
//    std::vector<TYPE> valarray;
//    for(int i=0;i<res*res;i++)
//    {
//        if(std::isnormal(values[i].real()))
//        {
//            valarray.push_back(values[i].real());
//        }
//    }
//    if(valarray.size()>0)
//    {
//        std::sort(valarray.begin(),valarray.end());
//        return valarray[valarray.size()/2];
//    }
//    else
//    {
//        return 0.0;
//    }
//}

QColor getColor(cplx value,TYPE & maxmod,TYPE & minmod,ColorMode mode)
{
    QColor color;

    bool flagmod=std::isnan(value.real()) || std::isnan(maxmod) || std::isnan(minmod) || (maxmod-minmod)==0;
    bool flagarg=std::isnan(value.imag());

    if(mode==MODE_MOD_ONLY && !flagmod)
    {
        double x=(value.real()-minmod)/(maxmod-minmod);
        double a=0;
        double b=255;
        double value_mod=a+(b-a)*x;

        color=QColor::fromHsl(0,0,clamp(value_mod,0,255));
    }
    else if(mode==MODE_MOD_ISO_ONLY && !flagmod)
    {
        double t_mod=log(1+value.real())*20;
        double x=std::pow(t_mod-int(t_mod),4);
        double a=128;
        double b=128-64;
        int value_mod=a+(b-a)*x;

        color=QColor::fromHsl(0,0,clamp(value_mod,0,255));
    }
    else if(mode==MODE_ARG_ONLY && !flagarg)
    {
        int value_arg=(value.imag()*180)/M_PI+180;

        color=QColor::fromHsl(clamp(value_arg,0,359),255,128);
    }
    else if(mode==MODE_ARG_AND_MOD_LIN && !flagarg && !flagmod)
    {
        double x=(value.real()-minmod)/(maxmod-minmod);
        double a=128-64;
        double b=128+64;
        int value_mod=a+(b-a)*x;

        int value_arg=(value.imag()*180)/M_PI+180;
        if(value_arg==360)value_arg=0;
        else if(value_arg>359)value_arg=359;
        else if(value_arg<0)value_arg=0;

        color=QColor::fromHsl(clamp(value_arg,0,359),255,clamp(value_mod,0,255));
    }
    else if(mode==MODE_ARG_AND_MOD_ISO && !flagarg && !flagmod)
    {
        double t_mod=log(1+value.real())*20;
        double x=std::pow(t_mod-int(t_mod),4);
        double a=128;
        double b=128-64;
        int value_mod=a+(b-a)*x;

        int value_arg=(value.imag()*180)/M_PI+180;
        if(value_arg==360)value_arg=0;
        else if(value_arg>359)value_arg=359;
        else if(value_arg<0)value_arg=0;

        color=QColor::fromHsl(clamp(value_arg,0,359),255,clamp(value_mod,0,255));
    }
    else
    {
        color=QColor::fromRgb( 0,0,0 );
    }

    return color;
}

void build_image(QString file,Resolution res,cplx * values,TYPE & maxmod,TYPE & minmod,ColorMode mode)
{
    QImage image(res.x,res.y,QImage::Format_RGB32);
    for(unsigned int i=0;i<res.x;i++)
    {
        for(unsigned int j=0;j<res.y;j++)
        {
            int id=i+res.x*j;
            image.setPixel(i,res.y-1-j,getColor(values[id],maxmod,minmod,mode).rgb());
        }
    }
    image.save(file);
}

void paste_at(int x,int y,QImage & ia,QImage & ib)
{
    for(int i=0;i<ib.width();i++)
    {
        for(int j=0;j<ib.height();j++)
        {
            ia.setPixel(x+i,y+j,ib.pixel(i,j));
        }
    }
}

void plot_compile(QString filename,QString dir,int nb_per_col)
{
    QStringList list=QDir(dir).entryList(QDir::Files|QDir::NoDotAndDotDot);

    std::cout<<list.size()<<std::endl;
    QImage img_compile;

    for(int i=0;i<list.size();i++)
    {
        QImage img(dir+list[i]);

        if(i==0)
        {
            img_compile=QImage((list.size()/nb_per_col)*img.width(),nb_per_col*img.height(),QImage::Format_RGB32);
        }
        std::cout<<i<<" "<<img.width()<<" "<<img.height()<<std::endl;

        int x=(i/nb_per_col)*img.width();
        int y=(i%nb_per_col)*img.height();
        paste_at(x,y,img_compile,img);
    }

    img_compile.save(filename);
}
