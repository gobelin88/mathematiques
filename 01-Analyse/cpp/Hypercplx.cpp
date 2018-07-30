#include "Hypercplx.h"

void bench_hypercplx()
{
    std::cout<<"------------------------------------------------------"<<std::endl;
    HyperCplx<-1> a(1,2);
    HyperCplx<-1> b(2,1);
    std::cout<<a<<" "<<arg(a)<<" "<<abs(a)<<std::endl;
    std::cout<<a*b<<" "<<b*a<<" "<<b+a<<" "<<a+b<<std::endl;
    std::cout<<a*2.0<<" "<<2.0*a<<" "<<b+2.0*a<<" "<<a*2.0+b<<std::endl;
    std::cout<<a+2.0<<" "<<2.0+a<<" "<<b+2.0+a<<" "<<a+2.0+b<<std::endl;
    std::cout<<a-2.0<<" "<<2.0-a<<" "<<b-2.0-a<<" "<<a-2.0-b<<std::endl;
    std::cout<<-a<<std::endl;
    std::cout<<a/b<<" "<<a/2.0<<" "<<2.0/a<<std::endl;

    std::cout<<"------------------------------------------------------"<<std::endl;
    std::complex<double> A(1,2);
    std::complex<double> B(2,1);
    std::cout<<A<<" "<<arg(A)<<" "<<abs(A)<<std::endl;
    std::cout<<A*B<<" "<<B*A<<" "<<B+A<<" "<<A+B<<std::endl;
    std::cout<<A*2.0<<" "<<2.0*A<<" "<<B+2.0*A<<" "<<A*2.0+B<<std::endl;
    std::cout<<A+2.0<<" "<<2.0+A<<" "<<B+2.0+A<<" "<<A+2.0+B<<std::endl;
    std::cout<<A-2.0<<" "<<2.0-A<<" "<<B-2.0-A<<" "<<A-2.0-B<<std::endl;
    std::cout<<-A<<std::endl;
    std::cout<<A/B<<" "<<A/2.0<<" "<<2.0/A<<std::endl;
}
