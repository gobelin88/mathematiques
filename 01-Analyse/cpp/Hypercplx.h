#include <iostream>
#include <complex>

#ifndef HYPERCPLX_H
#define HYPERCPLX_H

template<int t>
class HyperCplx
{
public:
    HyperCplx()
    {
        this->a=0;
        this->b=0;
    }

    HyperCplx(double a)
    {
        this->a=a;
        this->b=0;
    }

    HyperCplx(double a,double b)
    {
        this->a=a;
        this->b=b;
    }

    double real()const{return a;}
    double imag()const{return b;}

    HyperCplx inv()const
    {
        double f=a*a-b*b*t;
        return HyperCplx(a/f,-b/f);
    }

    bool operator==(const HyperCplx & N)
    {
        return (a==N.real() && b==N.imag());
    }

    HyperCplx operator=(const HyperCplx & N)
    {
        this->a=N.a;
        this->b=N.b;
        return (*this);
    }

    HyperCplx operator+=(const HyperCplx & N)
    {
        this->a+=N.a;
        this->b+=N.b;
        return (*this);
    }

    HyperCplx operator*=(const HyperCplx & N)
    {
        *this=HyperCplx(this->real()*N.real()+t*this->imag()*N.imag(),
                        this->real()*N.imag()+this->imag()*N.real());
        return (*this);
    }

    HyperCplx operator-=(const HyperCplx & N)
    {
        this->a-=N.a;
        this->b-=N.b;
        return (*this);
    }

    HyperCplx operator-()
    {
        return HyperCplx<t>(-this->a,-this->b);
    }

    friend std::ostream& operator<<(std::ostream& os, const HyperCplx & N)
    {
        os<<"("<<N.real()<<","<<N.imag()<<")";
        return os;
    }

    //------------------------------
    HyperCplx exponential()const
    {
        if(t==-1)
        {
            double mod=exp(this->real());
            return HyperCplx(cos(this->imag())*mod,sin(this->imag())*mod);
        }
        else if(t==1)
        {
            double mod=exp(this->real());
            return HyperCplx(cosh(this->imag())*mod,sinh(this->imag())*mod);
        }
        else if(t==0)
        {
            double mod=exp(this->real());
            return HyperCplx(mod,mod*this->imag());
        }
    }

private:
    double a,b;
};

template<int t>
HyperCplx<t> operator*(const HyperCplx<t> & Na,const HyperCplx<t> & Nb)
{
    return HyperCplx<t>(Na.real()*Nb.real()+t*Na.imag()*Nb.imag(),
                        Na.real()*Nb.imag()+Na.imag()*Nb.real());
}

template<int t>
HyperCplx<t> operator*(const HyperCplx<t> & Na,const double & Nb)
{
    return HyperCplx<t>(Na.real()*Nb,Na.imag()*Nb);
}

template<int t>
HyperCplx<t> operator*(const double & Nb,const HyperCplx<t> & Na)
{
    return HyperCplx<t>(Na.real()*Nb,Na.imag()*Nb);
}

template<int t>
HyperCplx<t> operator+(const HyperCplx<t> & Na,const double & Nb)
{
    return HyperCplx<t>(Na.real()+Nb,Na.imag());
}

template<int t>
HyperCplx<t> operator/(const HyperCplx<t> & Na,const double & Nb)
{
    return HyperCplx<t>(Na.real()/Nb,Na.imag()/Nb);
}

template<int t>
HyperCplx<t> operator/(const double & Nb,const HyperCplx<t> & Na)
{
    return Nb*Na.inv();
}

template<int t>
HyperCplx<t> operator/(const HyperCplx<t> & Na,const HyperCplx<t> & Nb)
{
    return Na*Nb.inv();
}

template<int t>
HyperCplx<t> operator+(const double & Nb,const HyperCplx<t> & Na)
{
    return HyperCplx<t>(Na.real()+Nb,Na.imag());
}

template<int t>
HyperCplx<t> operator-(const HyperCplx<t> & Na,const double & Nb)
{
    return HyperCplx<t>(Na.real()-Nb,Na.imag());
}

template<int t>
HyperCplx<t> operator-(const double & Nb,const HyperCplx<t> & Na)
{
    return HyperCplx<t>(Nb-Na.real(),-Na.imag());
}

template<int t>
HyperCplx<t> operator+(const HyperCplx<t> & Na,const HyperCplx<t> & Nb)
{
    return HyperCplx<t>(Na.real()+Nb.real(),Na.imag()+Nb.imag());
}

template<int t>
HyperCplx<t> operator-(const HyperCplx<t> & Na,const HyperCplx<t> & Nb)
{
    return HyperCplx<t>(Na.real()-Nb.real(),Na.imag()-Nb.imag());
}

template<int t>
double arg(const HyperCplx<t> & N)
{
    return atan2(N.imag(),N.real());
}

template<int t>
double abs(const HyperCplx<t> & N)
{
    return sqrt(N.real()*N.real()+N.imag()*N.imag());
}

void bench_hypercplx();

#endif // HYPERCPLX_H
