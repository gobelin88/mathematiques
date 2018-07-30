#include "dl.h"

//--------------------------------------------------------------------
//Developpements Limités
//--------------------------------------------------------------------
cplx T_coef(cplx z,cplx a,const std::vector<cplx> & C)
{
    cplx sum=0.0,sump=0.0;
    unsigned int n=0;
    do
    {
        sump=sum;
        sum+=C[n]*func_pow(z-a,n);
        n++;
    }
    while(n<C.size());

    return sum;
}

//-----------------------------------------------------------------------
//f=1/( (1-z)(2-z) )

cplx func_poles_1_2(cplx z)
{
    return cplx(1,0)/( ( cplx(1,0)-z )*(cplx(2,0)-z ) );
}

cplx L_poles_1_2(cplx z,int it)
{
    cplx sum(0,0);
    for(int n=1;n<it;++n)
    {
        sum+=(func_pow(2.0,n-1.0)-1.0)/func_pow(z,n);
    }
    return sum;
}

cplx T_poles_1_2(cplx z,int it)
{
    cplx sum(0,0);
    for(int n=0;n<it;++n)
    {
        sum+=(1.0-1.0/func_pow(2,n+1))*func_pow(z,n);
    }
    return sum;
}

cplx TL_poles_1_2(cplx z,int it)
{
    cplx sum(0,0);
    for(int n=-it;n<it;++n)
    {
        if(n<0)
        {
            sum+=-TYPE(1.0)/func_pow(z,-n);
        }
        else
        {
            sum+= TYPE(-1.0)/func_pow(2,n+1) * func_pow(z,n);
        }
    }
    return sum;
}

//-----------------------------------------------------------------------
//f=1/(1-z)
cplx func_geometric(cplx z)
{
    cplx one=cplx(1,0);
    return one/(one-z);
}

cplx TorL_geometric(cplx z,int it)
{
    if(abs(z)<1.0)
    {
        return T_geometric(z,it);
    }
    else
    {
        return L_geometric(z,it);
    }
}

cplx L_geometric(cplx z,int it)   //1/(1-z)
{
    if(abs(z)>0)
    {
        cplx sum=1.0,sump=0.0;
        unsigned int n=1;
        do
        {
            sump=sum;
            sum+=func_pow(TYPE(1.0)/z,n);
            n++;
        }
        while((int)n<it);
        return -sum/z;
    }
    else
    {
        return cplx(0,0);
    }
}

cplx T_geometric(cplx z,int it)
{
    cplx sum=0.0,sump=0.0;
    unsigned int n=0;
    do
    {
        sump=sum;
        sum+=func_pow(z,n);
        n++;
    }
    while((int)n<it);
    return sum;
}
cplx T_geometric_dm(cplx z,int it,int m)
{
    cplx sum=0.0,sump=0.0;
    unsigned int n=m;
    do
    {
        sump=sum;
        sum+=func_pow(z,n-m)*nb_down_product(n,m);
        n++;
    }
    while((int)n<it);
    return sum;
}

std::vector<cplx> T_coef_geometric(int it,int it2,cplx a)
{
    std::vector<cplx> C(it);
    for(int i=0;i<it;i++)
    {
        C[i]=T_geometric_dm(a,it2,i)/nb_factorial(i);
    }
    return C;
}

cplx TT_geomertic_ac(cplx z,int it)
{
    static cplx a(0.5,0.5);
    static std::vector<cplx> c=T_coef_geometric(it,1000,a);
    return T_coef(z,a,c);
}

//--------------------------------------------------------------------
//f=1/(1+z)
cplx T_geometric_altern(cplx z,int it)   //1/(1+z)=
{
    cplx sum=0.0,sump=0.0;
    unsigned int n=0;

    do
    {
        sump=sum;
        sum+=func_pow(z,n)*nb_minus_one_power_n(n);
        n++;
    }
    while((int)n<it);
    return sum;
}

//--------------------------------------------------------------------
//f=ln(1+z)
cplx T_ln(cplx z,int it)   //ln(1+z)=
{
    cplx sum=0.0;
    unsigned int n=1;

    do
    {
        sum-=func_pow(z,n)/cplx(n*nb_minus_one_power_n(n),0);
        n++;
    }
    while((int)n<it);
    return sum;
}

//--------------------------------------------------------------------
//f=exp(z)
cplx T_exp(cplx z,int it)   //ln(1+z)=
{
    cplx sum=1.0;
    unsigned int n=1;

    do
    {
        sum+=func_pow(z,n)/nb_factorial(n);
        n++;
    }
    while((int)n<it);
    return sum;
}

//--------------------------------------------------------------------
//f=sin(z)
cplx T_sin(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=0;
    do
    {
        sum+=func_pow(z,1+2*n)/nb_factorial(1+2*n)*nb_minus_one_power_n(n);
        n++;
    }
    while((int)n<it);

    return sum;
}

//--------------------------------------------------------------------
//f=sinh(z)
cplx T_sinh(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=0;
    do
    {
        sum+=func_pow(z,1+2*n)/nb_factorial(1+2*n);
        n++;
    }
    while((int)n<it);

    return sum;
}

//--------------------------------------------------------------------
//f=cos(z)
cplx T_cos(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=0;
    do
    {
        sum+=func_pow(z,2*n)/nb_factorial(2*n)*nb_minus_one_power_n(n);
        n++;
    }
    while((int)n<it);

    return sum;
}

//--------------------------------------------------------------------
//f=cosh(z)
cplx T_cosh(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=0;
    do
    {
        sum+=func_pow(z,2*n)/nb_factorial(2*n);
        n++;
    }
    while((int)n<it);

    return sum;
}

//-----------------------------------------------------------------------
//f=z/(e^-1)
cplx func_bernoulli(cplx z)
{
    return z/(func_exp(z)-TYPE(1.0));
}

cplx T_bernoulli(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=0;
    do
    {
        sum+=(nb_bernoulli(n)*func_pow(z,n))/nb_factorial(n);
        n++;
    }
    while((int)n<it);

    return sum;
}

//--------------------------------------------------------------------
//f=tan(z)
cplx T_tan(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=1;
    do
    {
        sum+=((func_pow(2,2*n)-func_pow(2,4*n))*nb_bernoulli(2*n))/nb_factorial(2*n)*nb_minus_one_power_n(n)*func_pow(z,2*n-1.0);
        n++;
    }    
    while((int)n<it);

    return sum;
}

//--------------------------------------------------------------------
//f=tanh(z)
cplx T_tanh(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=1;
    do
    {
        sum+=(func_pow(2,4*n)-func_pow(2,2*n))*nb_bernoulli(2*n)/nb_factorial(2*n)*func_pow(z,2*n-1);
        n++;
    }
    while((int)n<it);

    return sum;
}

//--------------------------------------------------------------------
//f=cot(z)
cplx T_cot(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=0;
    do
    {
        sum+=func_pow(2,2*n)*nb_bernoulli(2*n)/nb_factorial(2*n)*func_pow(z,2*n)*nb_minus_one_power_n(n);
        n++;
    }
    while((int)n<it);

    return sum/z;
}

//--------------------------------------------------------------------
//f=coth(z)
cplx T_coth(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=0;
    do
    {
        sum+=func_pow(2,2*n)*nb_bernoulli(2*n)/nb_factorial(2*n)*func_pow(z,2*n);
        n++;
    }
    while((int)n<it);

    return sum/z;
}

//--------------------------------------------------------------------
//f=sqrt(z)
cplx T_sqrt(cplx z,int it)
{
    cplx sum=0.0;
    unsigned int n=0;
    do
    {
        sum+=nb_binomial(.5,n)*func_pow(z,n);
        n++;
    }
    while((int)n<it);

    return sum;
}
