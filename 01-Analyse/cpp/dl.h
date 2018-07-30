#include "typedef.h"
#include "nombres.h"
#include <vector>
#include <iostream>
#include <functions.h>

//--------------------------------------------------------------------
//Developpements Limités
//--------------------------------------------------------------------
cplx T_coef(cplx z,cplx a,const std::vector<cplx> & C);

//-----------------------------------------------------------------------
//f=1/( (1-z)(2-z) )
cplx func_poles_1_2(cplx z);
cplx L_poles_1_2(cplx z,int it);   //|z|>2
cplx T_poles_1_2(cplx z,int it);   //|z|<1
cplx TL_poles_1_2(cplx z,int it);  //1<|z|<2

//-----------------------------------------------------------------------
//f=1/(1-z)
cplx func_geometric(cplx z);
cplx TorL_geometric(cplx z,int it);
cplx L_geometric(cplx z,int it);
cplx T_geometric(cplx z,int it);
cplx T_geometric_dm(cplx z,int it,int m);                   //Dérivée m de la série de Taylor
std::vector<cplx> T_coef_geometric(int it,int it2,cplx a);  //Taylor du de la série de Taylor en a
cplx TT_geomertic_ac(cplx z,int it);                        //Continuation Analytique en 0.5 0.5

//--------------------------------------------------------------------
//f=1/(1+z)
cplx T_geometric_altern(cplx z,int it);

//--------------------------------------------------------------------
//f=ln(1+z)
cplx T_ln(cplx z,int it);

//--------------------------------------------------------------------
//f=exp(z)
cplx T_exp(cplx z,int it);

//--------------------------------------------------------------------
//f=sin(z)
cplx T_sin(cplx z,int it);

//--------------------------------------------------------------------
//f=cos(z)
cplx T_cos(cplx z,int it);

//--------------------------------------------------------------------
//f=sinh(z)
cplx T_sinh(cplx z,int it);

//--------------------------------------------------------------------
//f=cosh(z)
cplx T_cosh(cplx z,int it);

//-----------------------------------------------------------------------
//f=z/(e^-1)
cplx func_bernoulli(cplx z);
cplx T_bernoulli(cplx z,int it);

//--------------------------------------------------------------------
//f=tan(z)
cplx T_tan(cplx z,int it);

//--------------------------------------------------------------------
//f=tanh(z)
cplx T_tanh(cplx z,int it);

//--------------------------------------------------------------------
//f=cot(z)
cplx T_cot(cplx z,int it);

//--------------------------------------------------------------------
//f=coth(z)
cplx T_coth(cplx z,int it);

//--------------------------------------------------------------------
//f=sqrt(1+z)
cplx T_sqrt(cplx z,int it);

