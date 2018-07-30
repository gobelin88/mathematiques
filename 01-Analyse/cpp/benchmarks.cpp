#include "benchmarks.h"

void padding(std::string & str)
{
    if(str.size()<PADDDING_SIZE)
    {
        int delta=PADDDING_SIZE-(unsigned int)str.size();

        for(int i=0;i<delta;i++)
        {
            str.append(" ");
        }
    }
}

bool approx(double a,double b)
{
    if(abs(a+b)>0)
    {
        if(abs((a-b)/((a+b)*0.5))<PRECISION){return true;}
        else{return false;}
    }
    else
    {
        return true;
    }
}
bool approx(cplx a,cplx b)
{
    if(abs(a+b)>0)
    {
        if( abs((a-b)/((a+b)*0.5))<PRECISION ){return true;}
        else{return false;}
    }
    else
    {
        return true;
    }
}
void test(cplx a,cplx b,std::string str1,std::string str2)
{
    padding(str1);
    padding(str2);

    std::cout<<std::setprecision(PRECISION_N)<<std::setw(WIDTH)<<std::setfill (' ')<<str1<<"=";
    std::cout<<std::setprecision(PRECISION_N)<<std::setw(WIDTH)<<std::setfill (' ')<<a;
    std::cout<<std::setprecision(PRECISION_N)<<std::setw(WIDTH)<<std::setfill (' ')<<str2<<"=";
    std::cout<<std::setprecision(PRECISION_N)<<std::setw(WIDTH)<<std::setfill (' ')<<b;

    if(approx(a,b))
        {std::cout<<"       (ok)"<<std::endl;}
    else{std::cout<<"       (failed) delta="<<log(abs((b-a)/a))/log(10)<<std::endl;}
}

void test(double a,double b,std::string str1,std::string str2)
{
    padding(str1);
    padding(str2);

    std::cout<<std::setprecision(PRECISION_N)<<std::setw(WIDTH)<<std::setfill (' ')<<str1<<"=";
    std::cout<<std::setprecision(PRECISION_N)<<std::setw(WIDTH)<<std::setfill (' ')<<a;
    std::cout<<std::setprecision(PRECISION_N)<<std::setw(WIDTH)<<std::setfill (' ')<<str2<<"=";
    std::cout<<std::setprecision(PRECISION_N)<<std::setw(WIDTH)<<std::setfill (' ')<<b;
    if(approx(a,b))
        {std::cout<<"       (ok)"<<std::endl;}
    else{std::cout<<"       (failed) delta="<<log(std::abs((b-a)/a))/log(10)<<std::endl;}
}


void bench_functions()
{
    std::cout<<"Bench functions-----------------------------------------------------------"<<std::endl;
    test(f_gamma(5,5000000)    ,f_gamma(4,5000000)*TYPE(4.0)   ,"f_gamma(5)      ","f_gamma(4)*4        ");
    test(f_gamma(0.5,100000)   ,sqrt(M_PI)               ,"f_gamma(0.5)    ","sqrt(pi)            ");
    test(f_zeta_form0(2,200000),M_PI*M_PI/6.0            ,"f_zeta_form0(2) ","pi^2/6              ");
    test(f_zeta_form1(2,9000)  ,M_PI*M_PI/6.0            ,"f_zeta_form1(2) ","pi^2/6              ");
    test(f_zeta_form2(2,18000) ,M_PI*M_PI/6.0            ,"f_zeta_form2(2) ","pi^2/6              ");
    test(f_zeta_form3(2,18000) ,M_PI*M_PI/6.0            ,"f_zeta_form3(2) ","pi^2/6              ");
    test(f_zeta_form4(2,18000) ,M_PI*M_PI/6.0            ,"f_zeta_form4(2) ","pi^2/6              ");
    test(f_zeta_form5(2,4500)  ,M_PI*M_PI/6.0            ,"f_zeta_form5(2) ","pi^2/6              ");
    test(f_theta(2,64)         ,sqrt(TYPE(0.5))*f_theta(0.5,64),"f_theta(2)      ","sqrt(1/2) theta(1/2)");
    test(f_xi_form0(0,1500)    ,0.5                       ,"f_xi_form0(0)   ","0.5                 ");
    test(f_xi_form0(1,1500)    ,0.5                       ,"f_xi_form0(1)   ","0.5                 ");
    test(f_xi_form0(2,1500)    ,M_PI/6.0                  ,"f_xi_form0(2)   ","-M_PI/6.0           ");
    test(f_xi_form1(0,1500)    ,0.5                       ,"f_xi_form1(0)   ","0.5                 ");
    test(f_xi_form1(1,1500)    ,0.5                       ,"f_xi_form1(1)   ","0.5                 ");
    test(f_xi_form1(2,1500)    ,M_PI/6.0                  ,"f_xi_form1(2)   ","-M_PI/6.0           ");
    std::cout<<std::endl;
    test(f_li(2)                 ,1.0451637801174927848445888891946131365226155781512,"f_li(2)        ",".        ");
    test(f_li(10)                ,6.1655995047872979375229817526695227491306028063765,"f_li(10)       ",".        ");
    test(f_li(30)                ,13.022632254541379652933606259562175628524320349270,"f_li(30)       ",".        ");
    test(f_li(50)                ,18.468696364806183438727930857327314052824834270989,"f_li(50)       ",".        ");
    test(f_li(100)               ,30.126141584079629925901741339032184979599907038905,"f_li(100)      ",".        ");
    test(f_li(200)               ,50.192171165963783809779123918235093014416305769480,"f_li(200)      ",".        ");
    test(f_li(500)               ,101.79387248862665948799642736109395432440039920759,"f_li(500)      ",".        ");
    test(f_li(cplx(1,1))     ,cplx(0.613911669221195508573360572483672790082596805739,
                                           2.059584214192577631503059203576739146551448640665),"f_li(1+i)      ",".        ");
    test(f_li(cplx(0,1))     ,cplx(0.472000651439568650777606107614127836507330543018,
                                           2.941558494949385099300999980021326772089446035251),"f_li(i)        ",".        ");
    test(f_li(cplx(0,-1))    ,cplx(0.472000651439568650777606107614127836507330543018,
                                           -2.941558494949385099300999980021326772089446035251),"f_li(-i)      ",".        ");
    std::cout<<std::endl;
    test(f_ei(1)                 ,1.8951178163559367554665209343316342690170605817327,"f_ei(1)        ",".        ");
    test(f_ei(10)                ,2492.2289762418777591384401439985248489896471014309,"f_ei(10)       ",".        ");
    test(f_ei(30)                ,3.6897320940727419706400632891084574699683612689e11,"f_ei(30)       ",".        ");
    test(f_ei(35)                ,4.6690550144661595445001462909900637974804574746e13,"f_ei(35)       ",".        ");
    test(f_ei(37)                ,3.2579889986722639967900016836981349125068665226e14,"f_ei(37)       ",".        ");
    test(f_ei(40)                ,6.0397182636112415783592314185106912937028885852e15,"f_ei(40)       ",".        ");
    test(f_ei(45)                ,7.9439160357044537715101683032183320016608854950e17,"f_ei(45)       ",".        ");
    test(f_ei(50)                ,1.0585636897131690963061541433229987195098919751e20,"f_ei(50)       ",".        ");
    test(f_ei(100)               ,2.7155527448538798219140146423108254102957939316e41,"f_ei(100)      ",".        ");
    test(f_ei(200)               ,3.6312352331593568523967100438464250464613074668e84,"f_ei(200)      ",".        ");
    test(f_ei(500)               ,2.812821397886294337474931517896438697693486786e214,"f_ei(500)      ",".        ");
    test(f_ei(-1)                ,-0.219383934395520273677163775460121649031047293406,"f_ei(-1)       ",".        ");
    test(f_ei(-10)               ,-4.15696892968532427740285981027818038434629008e-6,"f_ei(-10)       ",".        ");
    test(f_ei(-30)               ,-3.02155201068881254481582504515369792116735755e-15,"f_ei(-30)      ",".        ");
    test(f_ei(-35)               ,-1.75270593899473720005483266909291366670975558e-17,"f_ei(-35)      ",".        ");
    test(f_ei(-37)               ,-2.24702069758857122208335330790406424900840142e-18,"f_ei(-37)      ",".        ");
    test(f_ei(-40)               ,-1.03677326145165697215064188914525977128102461e-19,"f_ei(-40)      ",".        ");
    test(f_ei(-45)               ,-6.22569080946238364309530845334815389422949237e-22,"f_ei(-45)      ",".        ");
    test(f_ei(-50)               ,-3.78326402955045901869896785402128578030289318e-24,"f_ei(-50)      ",".        ");
    test(f_ei(-100)              ,-3.68359776168203218023519262050811898765522013e-46,"f_ei(-100)     ",".        ");
    test(f_ei(-200)              ,-6.88522610630763559771081748245579297383680869e-90,"f_ei(-200)     ",".        ");
    test(f_ei(-500)              ,-1.42207678225363842209819393605727828160786484e-220,"f_ei(-500)    ",".        ");

    //test(f_ei(1000)              ,1.972045137141238302809645048412023552690317566e431,"f_ei(1000)      ","0        ");
    test(f_ei(cplx(1,1))     ,cplx(1.764625985563854068426738161351237966008304411668,
                                           2.387769851510522419262792089103796064407333845441),"f_ei(1+i)      ",".        ");
    test(f_ei(cplx(0,1))     ,cplx(0.337403922900968134662646203889150769997578032585,
                                           2.516879397162079634172675005462931099910922654425),"f_ei(i)        ",".        ");
    test(f_ei(cplx(0,-1))    ,cplx(0.337403922900968134662646203889150769997578032585,
                                           -2.516879397162079634172675005462931099910922654425),"f_ei(-i)      ",".        ");
    std::cout<<std::endl;
    test(f_pi(100)       ,25    ,"f_pi(100)      ","25        ");
    std::cout<<"-------------------------------------------------------------------------------"<<std::endl;
}

void bench_TL_functions()
{
    std::cout<<"Bench_T_functions-------------------------------------------------------------"<<std::endl;

    test(T_geometric(0.5,1e3)  ,func_geometric(0.5),"T_geometric(0.5) ","func_geometric(0.5)");
    test(T_ln(0.5,1e3)         ,std::log(1+0.5)    ,"T_ln(0.5)        ","log(1+0.5)         ");
    test(T_exp(0.5,1e3)        ,std::exp(0.5)      ,"T_exp(0.5)       ","exp(0.5)           ");
    test(T_sinh(0.5,1e3)       ,std::sinh(0.5)     ,"T_sinh(0.5)      ","sinh(0.5)          ");
    test(T_cosh(0.5,1e3)       ,std::cosh(0.5)     ,"T_cosh(0.5)      ","cosh(0.5)          ");
    test(T_tanh(0.5,50)        ,std::tanh(0.5)     ,"T_tanh(0.5)      ","tanh(0.5)          ");
    test(T_sin(0.5,1e3)        ,std::sin(0.5)      ,"T_sin(0.5)       ","sin(0.5)           ");
    test(T_cos(0.5,1e3)        ,std::cos(0.5)      ,"T_cos(0.5)       ","cos(0.5)           ");
    test(T_tan(0.5,50)         ,std::tan(0.5)      ,"T_tan(0.5)       ","tan(0.5)           ");
    test(T_bernoulli(0.5,50)   ,func_bernoulli(0.5),"T_bernoulli(0.5) ","func_bernoulli(0.5)");
    test(T_poles_1_2(0.5,100)  ,func_poles_1_2(0.5),"T_poles_1_2(0.5) ","func_poles_1_2(0.5)");
    test(L_poles_1_2(2.5,100)  ,func_poles_1_2(2.5),"L_poles_1_2(2.5) ","func_poles_1_2(2.5)");
    test(TL_poles_1_2(1.5,100) ,func_poles_1_2(1.5),"TL_poles_1_2(1.5)","func_poles_1_2(1.5)");
    test(T_sqrt(0.5,10)        ,sqrt(0.5+1)        ,"T_sqrt(0.5)      ","sqrt(1.5)         ");
    test(T_sqrt(0.1,10)        ,sqrt(0.1+1)        ,"T_sqrt(0.1)      ","sqrt(1.1)         ");

    std::cout<<"-------------------------------------------------------------------------------"<<std::endl;
}

void bench_numbers()
{
    std::cout<<"Bench_numbers ------------------------------------------------------------------"<<std::endl;
    test(nb_factorial(5)         ,120 ,"nb_factorial(5)        ","120        ");
    test(nb_factorial(100)       ,9.332621544394415268169923885626670049071596826438e157 ,"nb_factorial(100)      ","120        ");
    test(nb_factorial(150)       ,5.713383956445854590478932865261054003189553578601e262 ,"nb_factorial(200)      ","120        ");
    test(nb_bernoulli(0)         ,1        ,"nb_bernoulli(0)        ","1            ");
    test(nb_bernoulli(1)         ,-1.0/2.0 ,"nb_bernoulli(1)        ","-1/2         ");
    test(nb_bernoulli(2)         ,1.0/6.0  ,"nb_bernoulli(2)        ","1/6          ");
    test(nb_bernoulli(3)         ,0        ,"nb_bernoulli(3)        ","0            ");
    test(nb_bernoulli(4)         ,-1.0/30.0,"nb_bernoulli(4)        ","-1/30        ");
    test(nb_minus_one_power_n(0) ,1        ,"nb_minus_one_power_n(0)","1            ");
    test(nb_minus_one_power_n(1) ,-1       ,"nb_minus_one_power_n(1)","-1           ");
    test(nb_binomial(5,3)        ,10       ,"nb_binomial(5,3)       ","10           ");
    test(nb_binomial(50,3)       ,19600    ,"nb_binomial(50,3)      ","19600        ");
    test(nb_mu(0)                ,0        ,"nb_mu(0)               ","0            ");
    test(nb_mu(1)                ,1        ,"nb_mu(1)               ","1            ");
    test(nb_mu(2)                ,-1       ,"nb_mu(2)               ","-1           ");
    test(nb_mu(3)                ,-1       ,"nb_mu(3)               ","-1           ");
    test(nb_mu(4)                ,0        ,"nb_mu(2)               ","0            ");
    test(nb_mu(5)                ,-1       ,"nb_mu(5)               ","-1           ");
    test(nb_mu(137)              ,-1       ,"nb_mu(137)             ","-1           ");

    std::cout<<"--------------------------------------------------------------------------------"<<std::endl;
}
