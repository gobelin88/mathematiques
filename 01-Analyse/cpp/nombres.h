#include <vector>

#include "typedef.h"
#include "primes.h"
#include "bernoulli.h"
#include "zeta.h"

#ifndef NOMBRES_H
#define NOMBRES_H

TYPE nb_down_product(int n,int m);                  //n*(n-1)*(n-2)*...*(n-m)   //(ok)
TYPE nb_factorial(TYPE n);                          //2*3*4*...*n               //(ok)
TYPE nb_minus_one_power_n(unsigned int n);          //(-1)^n                    //(ok)
TYPE nb_binomial(TYPE n,unsigned int k);            //C(n,k)=n!/(k!(n-k)!)      //(ok)
TYPE nb_bernoulli(unsigned int n);                  //Bn                        //(ok)

unsigned int nb_primes(unsigned int n);    //                          //(ok)
unsigned int nb_smallest_divisor(unsigned int n);
int nb_mu(unsigned int n);                 //Fonction de Moebius

#endif // NOMBRES_H

