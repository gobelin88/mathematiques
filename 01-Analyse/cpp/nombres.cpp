#include "nombres.h"

TYPE nb_down_product(int n,int m)    //n*(n-1)*(n-2)*...*(n-m)
{
   if(m>0)
   {
    TYPE p=n;
    for(int k=1;k<m;k++){p*=(n-k);}
    return p;
   }
   else
   {
       return 1;
   }
}

TYPE nb_factorial(TYPE n)       //2*3*4*...*n
{
    if (n<200)
    {
        TYPE value=1;
        for(unsigned int i=1;i<=n;i+=1){value*=i;}
        return value;
    }
    else
    {
        return TYPE_MAX;
    }
}

TYPE nb_minus_one_power_n(unsigned int n)
{
    if(n&1){return TYPE(-1);}
    else {return TYPE(1);}
}

TYPE nb_binomial(TYPE n,unsigned int k)
{
    //return nb_factorial(n)/(nb_factorial(k)*nb_factorial(n-k));
    TYPE prod=1;
    for(unsigned int j=1;j<=k;j++)
    {
        prod*=( TYPE(n+1)/TYPE(j)-1 );
    }
    return prod;
}

TYPE nb_bernoulli(unsigned int n)
{
    if(n<NB_BERNOULLI)
    {
        return bernoulli_list[n];
    }
    else if( (n&1)==1){return 0.0;}

    double sum=0.0;

    for(unsigned int k=0;k<n;k++)
    {
        sum+=nb_bernoulli(k)*nb_binomial(n+1,k);
    }

    return -1.0/(n+1.0)*sum;
}

unsigned int nb_primes(unsigned int n)
{
    return primes_list[n];
}

unsigned int nb_smallest_divisor(unsigned int n)
{
    unsigned int srx=(unsigned int)(sqrt(n)+1);
    if(n==1)
    {
        return 1;
    }
    else
    {
        for(int i=0;i<NB_PRIMES;i++)
        {
            if( n % primes_list[i]==0 ){return primes_list[i];}
            else if(primes_list[i]>srx){return n; }
            else{return 0;}
        }
        return 0;
    }
}

int nb_mu(unsigned int n)
{
    if(n==1)return 1;

    std::vector<unsigned int> dlist;
    unsigned int sr=n;
    do
    {
        unsigned int d=nb_smallest_divisor(sr);

        bool already=false;
        for(int i=0;i<dlist.size();i++)
        {
            if(dlist[i]==d)
            {
                already=true;
                break;
            }
        }
        if(already==false)dlist.push_back(d);
        else
        {
            return 0;
        }
        sr/=d;
    }
    while(sr!=1);

    return (int)nb_minus_one_power_n( (unsigned int) dlist.size() );
}
