#include "functions.h"

///////////////////////////////////////////////////////////////////////////////////////////
/// STD

double func_Arg(cplx z,double alpha)
{
    double a=std::atan2(z.imag(),z.real());
    double b=(a<0)?2*M_PI+a:a;


    if(b>=0 && b<alpha)
    {
        return b;
    }
    else
    {
        b=-(2*M_PI-b);

        return b;
    }
}
cplx func_Ln(cplx z){return func_Lna(z,M_PI);}
cplx func_Lna(cplx z,double alpha){return cplx(log(abs(z)),func_Arg(z,alpha));}
cplx func_Lnak(cplx z,double alpha,int k){return cplx(log(abs(z)),func_Arg(z,alpha)+k*2*M_PI);}
cplx func_Lnp(cplx z){return func_Ln(z+cplx(1,0));}

cplx func_exp(cplx z)
{
    return cplx(cos(z.imag()),sin(z.imag()))*exp(z.real());
    //return z.exponential();
}
cplx func_cos(cplx z)
{
    return (func_exp(z*cplx(0,1))+func_exp(-z*cplx(0,1)))/2.0;
}
cplx func_sin(cplx z)
{
    return (func_exp(z*cplx(0,1))-func_exp(-z*cplx(0,1)))/(2.0*cplx(0,1));
}
cplx func_cosh(cplx z)
{
    return (func_exp(z)+func_exp(-z))/2.0;
}
cplx func_sinh(cplx z)
{
    return (func_exp(z)-func_exp(-z))/2.0;
}

cplx func_tan(cplx z){return func_sin(z)/func_cos(z);}
cplx func_tanh(cplx z){return func_sinh(z)/func_cosh(z);}
cplx func_cot(cplx z){return TYPE(1.0)/func_tan(z);}
cplx func_coth(cplx z){return TYPE(1.0)/func_tanh(z);}

cplx func_pow(cplx z,cplx p)
{
    return func_exp(p*func_Ln(z));
}

cplx func_pow(cplx z,double p)
{
    return func_exp(p*func_Ln(z));
}

cplx func_pow(double z,cplx p)
{
    return func_exp(p*func_Ln(z));
}

cplx func_pow(double z,double p)
{
    return func_exp(p*func_Ln(z));
}

cplx func_sqrt(cplx z)
{
    return func_pow(z,0.5);
}

cplx f_gamma(cplx s,int it)
{
    cplx sum=func_pow(it,s)/s;
    int k=1;
    while(k<it)
    {
        sum*=TYPE(k)/(TYPE(k)+s);
        k++;
    }
    return sum;
}

cplx f_digamma(cplx s,int it)
{
    cplx sum=-EULAR_MACHERONIS;
    int k=1;
    while(k<it)
    {
        sum+=(s-TYPE(1.0))/(TYPE(k)*(TYPE(k)+s-TYPE(1.0)));
        k++;
    }
    return sum;
}

cplx f_beta(cplx x,cplx y,int it)
{
    return f_gamma(x,it)*f_gamma(y,it)/f_gamma(x+y,it);
}

cplx f_rbeta(cplx s,int it)
{
    return f_beta(s.real(),s.imag(),it);
}

cplx    f_wallis(cplx s,int it )
{
    return sqrt(M_PI)*TYPE(0.5)*f_gamma(s*TYPE(0.5)+TYPE(0.5),it)/f_gamma(s*TYPE(0.5)+TYPE(1.0),it);
}

cplx f_theta(cplx s,int it)
{
    return f_thetaq(func_exp(-s*M_PI),it);
}

cplx f_thetaq(cplx s,int it)
{
    cplx sum(0,0);

    if(abs(s)<1)
    {
        for(int n=-it;n<=it;n++)
        {
            sum+=func_pow( s , n*n );
        }
    }
    else
    {
        return cplx(INFINITY,INFINITY);
    }

    return sum;
}

cplx f_zeta_form0(cplx s,int it)
{
    if(s.real()<1)return cplx(INFINITY,INFINITY);
    cplx sum(1,0);
    int n=2;

    while(n<(it+2))
    {
        sum+=func_pow(n,-s);
        n++;
    }

    return sum;
}

cplx f_zeta_form1(cplx s,int it)
{
    if(s.real()<1)return cplx(INFINITY,INFINITY);
    cplx one=cplx(1,0);
    cplx sum(1,0);
    unsigned int n=0;

    while((int)n<it)
    {
        sum*=one/(one-func_pow(double(primes_list[n]),-s));
        n++;
    }

    return sum;
}

//it=1  1
//it=2  1-2^(-s)
//it=3  1-2^(-s)+3^(-s)...
cplx f_zeta_form2(cplx s,int it)
{
    if(s.real()<0)return cplx(INFINITY,INFINITY);
    cplx one=cplx(1,0);
    return f_eta(s,it)/(one-func_pow(2.0,one-s));
}

cplx f_zeta_form3(cplx s,int it)
{
    if(s==cplx(0.0,0.0))
    {
        return -0.5;
    }
    else if(s.real()<0.5)
    {
        return f_zeta_form2(TYPE(1.0)-s,it)*f_gamma((TYPE(1.0)-s)/TYPE(2.0),it)/f_gamma(s/TYPE(2.0),it)*func_pow(M_PI,s-TYPE(0.5));
    }
    else
    {
        return f_zeta_form2(s,it);
    }
}

cplx f_zeta_form4(cplx s,int it)
{
    if(s==cplx(0.0,0.0))
    {
        return -0.5;
    }
    else if(s.real()<0.5)
    {
        return f_zeta_form2(TYPE(1.0)-s,it)*func_sin(TYPE(0.5)*s*M_PI)*f_gamma(TYPE(1.0)-s,it)*func_pow(M_PI*2.0,s)/M_PI;
    }
    else
    {
        return f_zeta_form2(s,it);
    }
}

cplx f_zeta_form5(cplx s,int it)
{
    if(it>NB_NT_ZEROS)return cplx(0,0);
    return f_xi_form1(s,it)*func_pow(M_PI,s*TYPE(0.5))/((s-TYPE(1.0))*f_gamma(s*TYPE(0.5)+TYPE(1.0),it));
}

cplx    f_zeta_prime_form0(cplx s,int it)
{
    if(s.real()<1)return cplx(INFINITY,INFINITY);
    cplx sum(1,0);
    int n=1;

    while(n<(it+2))
    {
        sum+=-func_pow(n,-s)*log(TYPE(n));
        n++;
    }

    return sum;
}

cplx    f_zeta_prime_form1(cplx s,int it)
{
    if(s.real()<0)return cplx(INFINITY,INFINITY);
    cplx sum(1,0);
    int n=1;

    cplx A=-(log(TYPE(2.0))*func_pow(TYPE(2.0),TYPE(1.0)-s))/func_pow(TYPE(1.0)-func_pow(TYPE(2.0),TYPE(1.0)-s),TYPE(2.0));
    cplx B=TYPE(1.0)/(TYPE(1.0)-func_pow(TYPE(2.0),TYPE(1.0)-s));

    while(n<(it+2))
    {
        sum+= (A-B*func_Ln(TYPE(n)))*nb_minus_one_power_n(n+1)*func_pow(n,-s);
        n++;
    }

    return sum;
}

cplx f_eta(cplx s,int it)
{
    if(s.real()<0)return cplx(INFINITY,INFINITY);
    cplx sum(0,0);
    int n=1;
    while(n<(it+1))
    {
        sum-=nb_minus_one_power_n(n)*func_pow(n,-s);
        n++;
    }

    return sum;
}

double f_zeta_crit(double x,int it)
{
    return abs(f_zeta_form2(cplx(0.5,x),it));
}

cplx f_zeta_spec(cplx z,int it)
{
    cplx sum(0,0);
    for(int i=0;i<it;i++)
    {
        sum+=func_cos( zeta_nt_zero[i]*func_Ln(z) );
    }
    return sum;
}

cplx f_xi_hat(cplx s,int it)
{
    if(std::abs(s)==INFINITY || std::isnan(std::abs(s)))
    {
        return cplx(0,0);
    }
    else if(s.real()>0.5)
    {
        return func_pow(M_PI,-s*TYPE(0.5))*f_gamma(s*TYPE(0.5),it)*f_zeta_form2(s,it);
    }
    else
    {
        return f_xi_hat(TYPE(1.0)-s,it);
    }
}

cplx f_xi_form0(cplx s,int it)
{
    if(s==cplx(1.0,0))
    {
        return TYPE(0.5);
    }
    else if(s.real()>=0.5)
    {
        return TYPE(0.5)*s*(s-TYPE(1.0))*func_pow(M_PI,-s*TYPE(0.5))*f_gamma(s*TYPE(0.5),it)*f_zeta_form2(s,it);
    }
    else
    {
        return f_xi_form0(TYPE(1.0)-s,it);
    }
}

cplx f_xi_form1(cplx s,int it)
{
    TYPE B=.5*log(4.0*M_PI)-0.5*EULAR_MACHERONIS-1;

    cplx prod(1,0);
    cplx S=s*(s-TYPE(1.0));

    for(int i=0;i<it;i++)
    {
        //double Rho=(0.25+zeta_nt_zero[i]*zeta_nt_zero[i]);//p(1-p)
        cplx Rho=cplx(0.5,zeta_nt_zero[i])*cplx(0.5,-zeta_nt_zero[i]);

        prod*=(TYPE(1.0)+S/Rho)*func_exp(s/Rho);
    }

    return TYPE(0.5)*func_exp(B*s)*prod;
}


double f_isprime(double x)
{
    if(x<3)return 1.0;

    unsigned int rx=(unsigned int)x;
    unsigned int sx=(unsigned int)sqrt(x);


    unsigned int i=0;
    unsigned int k=0;
    unsigned int nb=primes_list[0];
    while(nb<=sx)
    {
        if(rx%nb==0)return 0.0;

        i++;
        if(i<NB_PRIMES)
        {
            nb=primes_list[i];
        }
        else
        {
            k++;
            nb=MAX_PRIME+k;
        }
    }
    return 1;
}



double f_pi(double x)
{
    unsigned int sum=0;

    while( (double)primes_list[sum++]<=x );

    return (double)sum-1;
}



cplx continuedFractionOfE1(const cplx& z)
{
    cplx cf(0.,0.);
    int k(120);
    while ( k > 0 )
    {
        TYPE rk(k--);
        cf = rk/(1.+rk/(z+cf));
    }
    return z+cf;
}
cplx ascendingSeriesOfE1(const cplx& z)
{
    cplx rs(cplx(1.,0.)), tr(rs);
    int k(2);
    while ( abs(tr) > TYPE_EPSILON * abs(rs) )
    {
        TYPE rk(k++);
        tr *= (1/rk-1.)*(z/rk);
        rs += tr;
    }
    return z*rs;
}
cplx E1z(const cplx& z)
{
    cplx rs;
    TYPE zr(z.real()), za(abs(z));
    if ( za < 10. || ( zr < 0 && za < 20. ) )
    {
        rs = -EULAR_MACHERONIS - func_Ln(z) + ascendingSeriesOfE1(z);
        if ( zr <= 0 && z.imag() == 0 ) rs -= cplx(0.,M_PI);
    }
    else
    {
        rs = func_exp(0.0-z)/continuedFractionOfE1(z);
        //if ( zr <= 0 && z.imag() == 0 ) rs -= cplx(0.,M_PI);
    }

    return rs;
}

cplx f_ei(cplx s)
{
    return -E1z(-s);


    //    if(abs(s)==0)
    //    {
    //        return exp(s)/s;
    //    }
    //    else if(abs(s)<40)
    //    {
    //        cplx sum=0.0;
    //        cplx delta=0.0;
    //        unsigned int n=1;

    //        do
    //        {
    //            //delta=pow(s,n)/(n*nb_factorial(n));
    //            delta=TYPE(1.0)/n;for(unsigned int k=1;k<=n;k++){delta*=s/TYPE(k);}
    //            sum+=delta;
    //            n++;
    //        }
    //        while( abs(delta)>TYPE_EPSILON );

    //        return EULAR_MACHERONIS+log( cplx(abs(s.real()),s.imag()) ) +sum;
    //    }
    //    else
    //    {
    //        cplx sum=1.0,sump=0.0,sump_coef=0.0;

    //        cplx coef=exp(s)/s;

    //        cplx delta=0.0;
    //        unsigned int n=1;
    //        bool err=false;
    //        do
    //        {
    //            //delta=nb_factorial(n)/pow(s,n);
    //            delta=1;for(unsigned int k=1;k<=n;k++){delta*=TYPE(k)/s;}

    //            sump=sum+delta;
    //            sump_coef=sump*coef;
    //            if( std::isfinite(sump_coef.real()) && std::isfinite(sump_coef.imag()) )
    //            {
    //                sum=sump;
    //                n++;
    //            }
    //            else
    //            {
    //                err=true;break;
    //            }
    //        }
    //        while( abs(delta)>TYPE_EPSILON );

    ////        if(err)
    ////        {
    ////            std::cout<<n<<" "<<sum<<" "<<exp(s)/s<<" "<<std::isfinite(sum.real())<<" "<<std::isfinite(sum.imag())<<std::endl;
    ////        }

    //        return coef*sum;
    //    }
}

cplx f_li(cplx s)
{
    return f_ei(func_Ln(s));
}

double f_ei_d(double t)
{
    return f_ei(cplx(t,0)).real();
}

double f_li_d(double t)
{
    return f_li(cplx(t,0)).real();
}

//////////////////////////////////////////////////////////////////////////////////

void f_primes_list(const char * filename,int size)
{
    int n=0;

    std::ofstream os(filename);

    time_t t1,t2;
    time(&t1);

    unsigned int cpt=1;
    os<<"#ifndef PRIMES_H"<<std::endl;
    os<<"#define PRIMES_H"<<std::endl;

    os<<"#define NB_PRIMES "<<size<<std::endl;

    os<<"static const unsigned int primes_list["<< size <<"]={"<<std::endl;
    while(n<size)
    {
        cpt++;
        if(f_isprime(cpt))
        {
            os<<cpt<<",";
            n++;
            if(n%50==0)os<<"\n";

            std::cout<<double(n)/(size)*100.0<<" "<<cpt<<"\r";
            std::cout.flush();
        }
    }
    os<<"};"<<std::endl;
    os<<"#define MAX_PRIME "<<cpt<<std::endl;
    os<<"#endif"<<std::endl;
    os.close();

    time(&t2);
    std::cout<<std::endl;
    std::cout<<"dt="<<(t2-t1)<<"s"<<std::endl;
}

cplx f_sum(cplx(*f)(cplx),unsigned int a,unsigned int b)
{
    cplx sum(0,0);

    for(unsigned int k=a;k<b;k++)
    {
        sum+=f(k);
    }

    return sum;
}

cplx f_eular_maclaurin_integrate(cplx(*f)(cplx),
                                 cplx(*fn)(cplx,int),
                                 unsigned int a,
                                 unsigned int b,
                                 unsigned int it)
{
    cplx sum=f_sum(f,a,b)+(f(b)-f(a))*TYPE(0.5);
    for(unsigned int n=0;n<it;++n)
    {
        sum-=nb_bernoulli(2*n)/nb_factorial(2*n)*(fn(a,2*n-1)-fn(b,2*n-1));
    }
    return sum;
}

cplx    test0(cplx s)
{
    if(s.real()<0 && s.imag()<0)     {return cplx(-1,-1);}
    else if(s.real()<0 && s.imag()>0){return cplx(-1,1);}
    else if(s.real()>0 && s.imag()>0){return cplx(1,1);}
    else if(s.real()>0 && s.imag()<0){return cplx(1,-1);}
    else
    {
        return cplx(0,0);
    }
}

cplx    test0q(cplx s)
{
    return test0(-func_Ln(s));
}

cplx    test1(cplx s)
{
    if(s.real()<1 && s.imag()<0)return cplx(-1,-1);
    else if(s.real()<1 && s.imag()>0)return cplx(-1,1);
    else if(s.real()>1 && s.imag()>0)return cplx(1,1);
    else if(s.real()>1 && s.imag()<0)return cplx(1,-1);
    else
    {
        return cplx(0,0);
    }
}

cplx    test1q(cplx s)
{
    return test1(TYPE(1.0)-func_Ln(s));
}

cplx test2 (cplx s)
{
    TYPE t=2*M_PI*cos(2*M_PI*s.real())*cos(2*M_PI*s.imag());
    if(t>0)
    {
        return s;
    }
    else
    {
        return -s;
    }

}

cplx t_mobius(cplx z,cplx a, cplx b,cplx c,cplx d)
{
    return (a*z+b)/(c*z+d);
}

cplx t_mobiusn(cplx z,cplx xi1,cplx xi2,cplx k,double n)
{
    if(xi2.real()==INFINITY)
    {
        return std::pow(k,n)*(z-xi1)+xi1;
    }
    else if (xi1==xi2)
    {
        return z/(1.0+n*k*z);
    }
    else
    {
        cplx p=(z-xi1)/(z-xi2)*std::pow(k,n);
        return (xi2*p-xi1)/(p-1.0);
    }
}

cplx t_log(cplx z,cplx a, cplx b,cplx c,cplx d)
{
    return (a*func_Ln(z)+b)/(c*func_Ln(z)+d);
}
