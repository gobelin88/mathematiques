#include "nombres.h"
#include "iostream"
#include <fstream>
#include <time.h>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////////
/// STD
double func_Arg(cplx z,double alpha);
cplx func_Ln(cplx z);
cplx func_Lna(cplx z,double alpha);
cplx func_Lnp(cplx z);//Ln(1+z)
cplx func_Lnak(cplx z,double alpha,int k);

cplx func_exp(cplx z);
cplx func_cos(cplx z);
cplx func_sin(cplx z);
cplx func_cosh(cplx z);
cplx func_sinh(cplx z);
cplx func_tan(cplx z);
cplx func_tanh(cplx z);
cplx func_cot(cplx z);
cplx func_coth(cplx z);
cplx func_pow(cplx z,cplx p);
cplx func_pow(cplx z,double p);
cplx func_pow(double z,cplx p);
cplx func_pow(double z,double p);
cplx func_sqrt(cplx z);

////////////////////////////////////////////////////////////////////////////////////////////
///Transformations conformes
////////////////////////////////////////////////////////////////////////////////////////////
cplx t_mobius(cplx z,cplx a, cplx b,cplx c,cplx d);
cplx t_mobiusn(cplx z, cplx xi1,cplx xi2,cplx k,double n);
cplx t_log(cplx z,cplx a, cplx b,cplx c,cplx d);

////////////////////////////////////////////////////////////////////////////////////////////
///Théorie des nombres
////////////////////////////////////////////////////////////////////////////////////////////
//GAMMA
cplx    f_gamma(cplx s,int it);                     //Fonction gamma
cplx    f_digamma(cplx s,int it);                   //Fonction digamma
cplx    f_beta(cplx x,cplx y,int it);           //Fonction beta
cplx    f_rbeta(cplx s,int it);                     //Fonction beta(s.r,s.i)

//GAMMA BETA DERIVED
cplx    f_wallis(cplx s,int it );

//THETA
cplx    f_thetaq(cplx s,int it);                    //|s|<1 sum s^(n^2)
cplx    f_theta(cplx s,int it);                     //R>1 f_theta( e^(-pi s) )
//ZETA
cplx    f_zeta_form0(cplx s,int it);                //R(s)>1  définition de base par la somme de l'inverse des entiers     : sum 1/n^s
cplx    f_zeta_form1(cplx s,int it);                //R(s)>1  définition par le produit de gauss sur les nombres premiers  : prod (1/(1-p^(-s))
cplx    f_zeta_form2(cplx s,int it);                //R(s)>0  définition avec la fonction eta                              : f_eta(s)/(1-2^(1-s))
cplx    f_zeta_form3(cplx s,int it);                //s!=1 définition avec la relation fonctionelle symétrique             :
cplx    f_zeta_form4(cplx s,int it);                //s!=1 définition avec la relation fonctionelle asymétrique            :
cplx    f_zeta_form5(cplx s,int it);                //s!=1 définition par le produit de Weiertrass                         :

//ZETA PRIME
cplx    f_zeta_prime_form0(cplx s,int it);    //R(s)>1  définition de base : sum (-ln(n)^k)/n^s
cplx    f_zeta_prime_form1(cplx s,int it);    //R(s)>0

//Misc
cplx    f_eta(cplx s,int it);                       //R(s)>0 fonction zeta alternée                                        : sum (-1)^n/n^s
double      f_zeta_crit(double x,int it);                   //Module sur ma ligne critique (0.5+i*s)
cplx    f_zeta_spec(cplx z,int it);                 //
//XI
cplx    f_xi_hat(cplx s,int it);                    //Fonction xi avec pôle en 0 et 1                                      :pi^(-s/2) gamma(s/2) zeta(s)
cplx    f_xi_form0(cplx s,int it);                  //Fonction Xi de Riemann par l'équation fonctionnelle                  :
cplx    f_xi_form1(cplx s,int it);                  //Fonction Xi de Riemann par les zéros non triviaux                    :

//PRIMES
cplx continuedFractionOfE1(const cplx& z);
cplx ascendingSeriesOfE1(const cplx& z);
cplx E1z(const cplx& z);


cplx    f_ei(cplx s);                               //Exponentielle intégrale
double      f_ei_d(double t);                               //Exponentielle intégrale
cplx    f_li(cplx s);                               //Logarithme intégrale
double      f_li_d(double t);                               //Logarithme intégrale

template <unsigned int it, unsigned int nbzeros=0>
cplx      f_r(cplx s)                               //Fonction R de Riemann
{
    cplx summ=0.0;
    for(int n=1;n<it;n++)
    {
        cplx sum=f_li(func_pow(s,TYPE(1.0)/n))-LN2;//;
        for(int k=0;k<nbzeros;k++)
        {
            cplx zp(TYPE(0.5)/TYPE(n), zeta_nt_zero[k]/TYPE(n));
            sum-=2.0*f_ei( zp*func_Ln(s) ).real();

            //cplx zn(TYPE(0.5)/TYPE(n),-zeta_nt_zero[k]/TYPE(n));
            //sum-=f_ei( zp*std::log(s) )+f_ei( zn*std::log(s) );
        }
        summ+=TYPE(nb_mu(n))/TYPE(n) * sum;
    }
    return summ;
}

template <unsigned int it, unsigned int nbzeros=0>
double      f_r_d(double t)
{
    return f_r<it,nbzeros>(cplx(t,0)).real();
}



double      f_pi(double x);                                 //Nombre de nombre Premiers <x



double      f_isprime(double x);                            //x est t'il premier
void        f_primes_list(const char *filename,int size);   //Construit la liste des nombres premiers

////////////////////////////////////////////////////////////////////////////////////////////
///Sommes et intégrales
////////////////////////////////////////////////////////////////////////////////////////////
cplx    f_sum(cplx(*f)(cplx),
                  unsigned int a,
                  unsigned int b);                          //somme de a à b
cplx    f_eular_maclaurin_integrate(cplx(*f)(cplx),
                                        cplx(*fn)(cplx, int),
                                        unsigned int a,
                                        unsigned int b,
                                        unsigned int it);   //intégrale de Maclaurin de a à b (fn dérivée n-iéme de f)

////////////////////////////////////////////////////////////////////////////////////////////
///Misc
////////////////////////////////////////////////////////////////////////////////////////////
cplx    test0(cplx s);
cplx    test0q(cplx s);
cplx    test1(cplx s);
cplx    test1q(cplx s);
cplx    test2(cplx s);
