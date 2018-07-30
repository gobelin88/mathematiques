#include "nombres.h"
#include "plot.h"
#include "dl.h"

#include <iomanip>

#define WIDTH 40
#define PRECISION 1e-10
#define PRECISION_N 10
#define PADDDING_SIZE 30

bool approx(double a,double b);
bool approx(cplx a,cplx b);
void test(cplx a,cplx b,std::string str1,std::string str2);
void test(double a, double b, std::string str1, std::string str2);
void padding(std::string & str);


void bench_functions();
void bench_TL_functions();
void bench_numbers();
