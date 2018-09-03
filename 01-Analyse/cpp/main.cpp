#include <iostream>
#include <QElapsedTimer>


#include "typedef.h"
#include "benchmarks.h"
#include "fractals.h"

void plot_fractals(Resolution res)
{
    //    ----------------------------------------------------------------------------------------------
    //    Ensembles de Julia/MandelBrot pour f(z)=z^2
    double step=0.1;
    int id=0;
    double s=1.5;
    double ds=20*step;

    Resolution res_sub(128);
    for(double a=-s-ds*step;a<(s+step/2-ds*step);a+=step)
    {
        for(double b=-s;b<(s+step/2);b+=step)
        {
            plot_img(QString("./plots/fractals/mini_julia/mini_julia_%1_%2_%3.png").arg(id,4,10,QChar('0')).arg(a).arg(b),
                     [=](cplx z)->cplx{return frac_julia(z,256,cplx(a,b));},res_sub,-1.5,-1.5,1.5,1.5,MODE_ARG_AND_MOD_LIN);
            id++;
        }
    }
    plot_compile("./plots/fractals/julia_map.png","./plots/fractals/mini_julia/",2*s/step+1);

    plot_img("./plots/fractals/mandelbrot.png",[](cplx z)->cplx{return frac_mandelbrot(z,1024);},res,-1.8,-1.3,0.8,1.3,MODE_ARG_AND_MOD_LIN);
    plot_img(QString("./plots/fractals/julia_%1_%2.png").arg(-0.70).arg(0.30),[](cplx z)->cplx{return frac_julia(z,1024,cplx(-0.70,0.30));},res,-1.5,-1.5,1.5,1.5,MODE_ARG_AND_MOD_LIN);
    plot_img(QString("./plots/fractals/julia_%1_%2.png").arg(-0.83).arg(0.20),[](cplx z)->cplx{return frac_julia(z,1024,cplx(-0.83,0.20));},res,-1.5,-1.5,1.5,1.5,MODE_ARG_AND_MOD_LIN);
    plot_img(QString("./plots/fractals/julia_%1_%2.png").arg(-0.63).arg(0.70),[](cplx z)->cplx{return frac_julia(z,1024,cplx(-0.63,0.70));},res,-1.5,-1.5,1.5,1.5,MODE_ARG_AND_MOD_LIN);
    plot_img(QString("./plots/fractals/julia_%1_%2.png").arg(-0.23).arg(0.70),[](cplx z)->cplx{return frac_julia(z,1024,cplx(-0.23,0.70));},res,-1.5,-1.5,1.5,1.5,MODE_ARG_AND_MOD_LIN);
    plot_img(QString("./plots/fractals/julia_%1_%2.png").arg(-0.13).arg(-0.70),[](cplx z)->cplx{return frac_julia(z,1024,cplx(-0.13,-0.70));},res,-1.5,-1.5,1.5,1.5,MODE_ARG_AND_MOD_LIN);
    plot_img(QString("./plots/fractals/julia_%1_%2.png").arg(-0.13).arg(-0.80),[](cplx z)->cplx{return frac_julia(z,1024,cplx(-0.13,-0.80));},res,-1.5,-1.5,1.5,1.5,MODE_ARG_AND_MOD_LIN);

    //----------------------------------------------------------------------------------------------

    //plot_img(QString("./plots/fractals/mandelbrot_exp.png"),[](complexd z)->complexd{return frac_mandelbrot_set([](complexd s)->complexd{return exp(s);},z,1024,100);},4096 ,-10,-10,10,10,MODE_FRACTAL);
    //plot_img(QString("./plots/fractals/mandelbrot_sin.png"),[](complexd z)->complexd{return frac_mandelbrot_set([](complexd s)->complexd{return sin(s);},z,1024,100);},4096 ,-10,-10,10,10,MODE_FRACTAL);
    //plot_img(QString("./plots/fractals/mandelbrot_cos.png"),[](complexd z)->complexd{return frac_mandelbrot_set([](complexd s)->complexd{return cos(s);},z,1024,100);},4096 ,-10,-10,10,10,MODE_FRACTAL);

    //plot_img(QString("./plots/fractals/mandelbrot_zeta.png"),[](complexd z)->complexd{return frac_mandelbrot_set([](complexd s)->complexd{return f_zeta_ext4(s,128);},z,128,1000);},1024 ,-60,-40,20,40,MODE_FRACTAL);
    //plot_img_it("./plots/fractals/zeta_ext3_%1.png"           ,f_zeta_ext3    ,512,512,512,-30,-30,30,30);

}

void plot_mobius(Resolution res,ColorMode mode)
{
    double s=2;

    plot_img(QString("./plots/mobius/mobius1_elliptique_inverse/mobiusn_flux.png") ,[=](cplx z)->cplx{return t_mobiusn(z,1,-1,-1,1)-z;},res,-s,-s,s,s,mode);
    plot_img(QString("./plots/mobius/mobius1_elliptique_inverse/mobiusn.png") ,[=](cplx z)->cplx{return t_mobiusn(z,1,-1,-1,1);},res,-s,-s,s,s,mode);

    plot_img(QString("./plots/mobius/mobius2_elliptique_reflexion/mobiusn_flux.png") ,[=](cplx z)->cplx{return t_mobiusn(z,0.5,INFINITY,-1.0,1)-z;},res,-s,-s,s,s,mode);
    plot_img(QString("./plots/mobius/mobius2_elliptique_reflexion/mobiusn.png") ,[=](cplx z)->cplx{return    t_mobiusn(z,0.5,INFINITY,-1.0,1);},res,-s,-s,s,s,mode);

    for(int i=0;i<10;++i)
    {
        double t=(2.0*M_PI*i)/10.0;
        plot_riemann(QString("./plots/mobius/mobius1_elliptique_inverse/mobiusn_sphere_%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobiusn(z,1,-1,-1,1);},res,M_PI/4,t,mode);
        plot_riemann(QString("./plots/mobius/mobius1_elliptique_inverse/mobiusn_flux_sphere_%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobiusn(z,1,-1,-1,1)-z;},res,M_PI/4,t,mode);

        plot_riemann(QString("./plots/mobius/mobius2_elliptique_reflexion/mobiusn_sphere_%1.png").arg(i) ,[=](cplx z)->cplx{return    t_mobiusn(z,0.5,INFINITY,-1.0,1);},res,M_PI/4,t,mode);
        plot_riemann(QString("./plots/mobius/mobius2_elliptique_reflexion/mobiusn_flux_sphere_%1.png").arg(i) ,[=](cplx z)->cplx{return    t_mobiusn(z,0.5,INFINITY,-1.0,1)-z;},res,M_PI/4,t,mode);
    }

//        plot_img(QString("./plots/mobius/mobius2_elliptique_reflexion/mobiusn_flux_%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobiusn(z,0.5,INFINITY,-1.0,n)-z;},res,-s,-s,s,s,mode);
//        plot_img(QString("./plots/mobius/mobius2_elliptique_reflexion/mobiusn_%1.png").arg(i) ,[=](cplx z)->cplx{return    t_mobiusn(z,0.5,INFINITY,-1.0,n);},res,-s,-s,s,s,mode);
//        plot_riemann(QString("./plots/mobius/mobius2_elliptique_reflexion/mobiusn_sphere_%1.png").arg(i) ,[=](cplx z)->cplx{return    t_mobiusn(z,0.5,INFINITY,-1.0,n);},res,M_PI/4,M_PI/4,mode);

//        plot_img(QString("./plots/mobius/mobius3_elliptique_inversion_reflexion/mobiusn_flux_%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobiusn(z,(1.0-I*sqrt(3.0))/2.0,(1.0+I*sqrt(3.0))/2.0,(-1.0+I*sqrt(3.0))/2.0,n)-z;},res,-s,-s,s,s,mode);
//        plot_img(QString("./plots/mobius/mobius3_elliptique_inversion_reflexion/mobiusn_%1.png").arg(i) ,[=](cplx z)->cplx{return    t_mobiusn(z,(1.0-I*sqrt(3.0))/2.0,(1.0+I*sqrt(3.0))/2.0,(-1.0+I*sqrt(3.0))/2.0,n);},res,-s,-s,s,s,mode);
//        plot_riemann(QString("./plots/mobius/mobius3_elliptique_inversion_reflexion/mobiusn_sphere_%1.png").arg(i) ,[=](cplx z)->cplx{return    t_mobiusn(z,(1.0-I*sqrt(3.0))/2.0,(1.0+I*sqrt(3.0))/2.0,(-1.0+I*sqrt(3.0))/2.0,n);},res,M_PI/4,M_PI/4,mode);

//        plot_img(QString("./plots/mobius/mobius4_elliptique_reflexion_inversion/mobiusn_flux_%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobiusn(z,(1.0+I*sqrt(3.0))/2.0,(1.0-I*sqrt(3.0))/2.0,(-1.0+I*sqrt(3.0))/2.0,n)-z;},res,-s,-s,s,s,mode);
//        plot_img(QString("./plots/mobius/mobius4_elliptique_reflexion_inversion/mobiusn_%1.png").arg(i) ,[=](cplx z)->cplx{return    t_mobiusn(z,(1.0+I*sqrt(3.0))/2.0,(1.0-I*sqrt(3.0))/2.0,(-1.0+I*sqrt(3.0))/2.0,n);},res,-s,-s,s,s,mode);

//        plot_img(QString("./plots/mobius/mobius5_parabolique/mobiusn_flux_%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobiusn(z,0.0,0.0,1,n)-z;},res,-s,-s,s,s,mode);
//        plot_img(QString("./plots/mobius/mobius5_parabolique/mobiusn_%1.png").arg(i) ,[=](cplx z)->cplx{return    t_mobiusn(z,0.0,0.0,1,n);},res,-s,-s,s,s,mode);

//        plot_img(QString("./plots/mobius/mobius6_hyperbolique/mobiusn_flux_%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobiusn(z,(-1.0+sqrt(3.0))/2.0,(-1.0-sqrt(3.0))/2.0,7.0-4.0*sqrt(3.0),n)-z;},res,-s,-s,s,s,mode);
//        plot_img(QString("./plots/mobius/mobius6_hyperbolique/mobiusn_%1.png").arg(i) ,[=](cplx z)->cplx{return    t_mobiusn(z,(-1.0+sqrt(3.0))/2.0,(-1.0-sqrt(3.0))/2.0,7.0-4.0*sqrt(3.0),n);},res,-s,-s,s,s,mode);
}

void plot_all(Resolution res,ColorMode mode)
{
    double N=10;

//    //1 Similitude
//    for(int i=0;i<=N;++i)
//    {
//        double theta=(2.0*M_PI*i)/N;//0 à 2*pi
//        double t=(i-5)*2/N;         //-1 à 1
//        double k=(i+1)*0.2;            //1/5 à 11/5

//        plot_riemann(QString("./plots/spheres/mobius/1_similitude_rotation_sphere%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobius(z,exp(I*theta/2.0),0,0,exp(-I*theta/2.0));},res,M_PI/4,M_PI/4,mode);
//        plot_img(QString("./plots/spheres/mobius/1_similitude_rotation_plan%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobius(z,exp(I*theta/2.0),0,0,exp(-I*theta/2.0));},res,-2,-2,2,2,mode);
//        plot_riemann(QString("./plots/spheres/mobius/1_similitude_translation_sphere%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobius(z,1,t,0,1);},res,M_PI/4,M_PI/4,mode);
//        plot_img(QString("./plots/spheres/mobius/1_similitude_translation_plan%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobius(z,1,t,0,1);},res,-2,-2,2,2,mode);
//        plot_riemann(QString("./plots/spheres/mobius/1_similitude_homothetie_sphere%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobius(z,k,0,0,1/k);},res,M_PI/4,M_PI/4,mode);
//        plot_img(QString("./plots/spheres/mobius/1_similitude_homothetie_plan%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobius(z,k,0,0,1/k);},res,-2,-2,2,2,mode);
//        plot_riemann(QString("./plots/spheres/mobius/1_similitude_sphere%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobius(z,k*exp(I*theta/2.0),t/k*exp(-I*theta/2.0),0,1/k*exp(-I*theta/2.0));},res,M_PI/4,M_PI/4,mode);
//        plot_img(QString("./plots/spheres/mobius/1_similitude_plan%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobius(z,k*exp(I*theta/2.0),t/k*exp(-I*theta/2.0),0,1/k*exp(-I*theta/2.0));},res,-2,-2,2,2,mode);
//    }

//    //2 Transformation Inverse
//    plot_img("./plots/spheres/mobius/2_inverse_plan.png",[](cplx z)->cplx{return t_mobiusn(z,-1,1,-1,1);},res,  -2,-2,2,2,mode);
//    plot_img("./plots/spheres/mobius/2_inverse_flux_plan.png",[](cplx z)->cplx{return z-t_mobiusn(z,-1,1,-1,1);},res,  -2,-2,2,2,mode);
//    for(int i=0;i<=N;++i)
//    {
//        double t=(2.0*M_PI*i)/N;
//        plot_riemann(QString("./plots/spheres/mobius/2_inverse_sphere%1.png").arg(i) ,[](cplx z)->cplx{return t_mobiusn(z,-1,1,-1,1);},res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/mobius/2_inverse_flux_sphere%1.png").arg(i) ,[](cplx z)->cplx{return z-t_mobiusn(z,-1,1,-1,1);},res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/mobius/2_inverse_it_sphere%1.png").arg(i) ,[=](cplx z)->cplx{return t_mobiusn(z,-1,1,-1,(i+0.0001)/N*2);},res,M_PI/4,M_PI/4,mode);
//    }

//    return;

//    //Fonctions Spheres
//    for(int i=0;i<N;++i)
//    {
//        double t=(2.0*M_PI*i)/N;

//        plot_riemann(QString("./plots/spheres/func/f_xi_hat_sphere%1.png").arg(i) ,[](cplx z)->cplx{return f_xi_hat(z,200);},res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/f_gamma_sphere%1.png").arg(i) ,[](cplx z)->cplx{return f_gamma(z,200);},res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/f_digamma_sphere%1.png").arg(i) ,[](cplx z)->cplx{return f_digamma(z,200);},res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/zeta_form3_sphere%1.png").arg(i) ,[](cplx z)->cplx{return f_zeta_form3(z,200);},res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/func_exp_sphere%1.png").arg(i) ,func_exp,res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/func_cos_sphere%1.png").arg(i) ,func_cos,res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/func_sin_sphere%1.png").arg(i) ,func_sin,res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/func_tan_sphere%1.png").arg(i) ,func_tan,res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/func_cosh_sphere%1.png").arg(i) ,func_cosh,res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/func_sinh_sphere%1.png").arg(i) ,func_sinh,res,M_PI/4,t,mode);
//        plot_riemann(QString("./plots/spheres/func/func_tanh_sphere%1.png").arg(i) ,func_tanh,res,M_PI/4,t,mode);
//    }
//    return;

    //TESTS ---------------------------------------------------------------------
    plot_img("./plots/test_mod_only.png",[](cplx z)->cplx{return z;},res,  -10,-10,10,10,MODE_MOD_ONLY);
    plot_img("./plots/test_arg_only.png",[](cplx z)->cplx{return z;},res,  -10,-10,10,10,MODE_ARG_ONLY);
    plot_img("./plots/test_arg_and_mod_lin.png",[](cplx z)->cplx{return z;},res,  -10,-10,10,10,MODE_ARG_AND_MOD_LIN);
    plot_img("./plots/test_arg_and_mod_iso.png",[](cplx z)->cplx{return z;},res,  -10,-10,10,10,MODE_ARG_AND_MOD_ISO);

    //MOBIUS
    //plot_mobius(res,mode);

    //FUNCTIONS ---------------------------------------------------------------------
    plot_img("./plots/bernoulli.png",func_bernoulli ,res,-10,-10,10,10,mode);
    plot_img("./plots/geom.png"     ,func_geometric ,res,-2,-2,2,2,mode);
    plot_img("./plots/exp.png"      ,func_exp       ,res,-10,-10,10,10,mode);
    plot_img("./plots/Lnp.png"      ,func_Lnp        ,res,-2,-2,2,2,mode);
    plot_img("./plots/cos.png"      ,func_cos       ,res,-10,-10,10,10,mode);
    plot_img("./plots/sin.png"      ,func_sin       ,res,-10,-10,10,10,mode);
    plot_img("./plots/cosh.png"     ,func_cosh      ,res,-10,-10,10,10,mode);
    plot_img("./plots/sinh.png"     ,func_sinh      ,res,-10,-10,10,10,mode);
    plot_img("./plots/tan.png"      ,func_tan       ,res,-3,-3,3,3,mode);
    plot_img("./plots/tanh.png"     ,func_tanh      ,res,-3,-3,3,3,mode);
    plot_img("./plots/cot.png"      ,func_cot       ,res,-4,-4,4,4,mode);
    plot_img("./plots/coth.png"     ,func_coth      ,res,-4,-4,4,4,mode);
    plot_img("./plots/sqrt.png"     ,func_sqrt      ,res,-2,-2,2,2,mode);
    plot_img("./plots/ei.png"       ,f_ei           ,res,-20,-20,20,20,mode);
    plot_img("./plots/li.png"       ,f_li           ,res,-20,-20,20,20,mode);
    //plot_img("./plots/r10.png"      ,f_r<10>        ,res,-20,-20,20,20,mode);

    //Log déterminations
    for(int k=-1;k<=1;k++)
    {
        for(int a=0;a<=8;a++)
        {
            plot_img(QString("./plots/Ln_%1_%2.png").arg(k).arg(a)
                     ,[=](cplx z)->cplx{return func_Lnak(z,a*M_PI/4,k);}
            ,Resolution(512),-5,-5,5,5,MODE_ARG_AND_MOD_ISO);
        }
    }

    //SPECIALES-----------------------------------------------------------------
    plot_img_it("./plots/thetaq_%1.png"     ,f_thetaq       ,res,64,64,-1,-1,1,1,mode);
    plot_img_it("./plots/theta_%1.png"      ,f_theta        ,res,64,64,-1,-1,1,1,mode);



    plot_img_it("./plots/zeta_form0_%1.png" ,f_zeta_form0  ,res,512,512,-30,-30,30,30,mode);
    plot_img_it("./plots/zeta_form1_%1.png" ,f_zeta_form1  ,res,512,512,-30,-30,30,30,mode);
    plot_img_it("./plots/zeta_form2_%1.png" ,f_zeta_form2  ,res,512,512,-30,-30,30,30,mode);
    plot_img_it("./plots/zeta_form3_%1.png" ,f_zeta_form3  ,res,512,512,-30,-30,30,30,mode);
    plot_img_it("./plots/zeta_form4_%1.png" ,f_zeta_form4  ,res,512,512,-30,-30,30,30,mode);
    plot_img_it("./plots/zeta_form5_%1.png" ,f_zeta_form5  ,res,512,512,-30,-30,30,30,mode);
    plot_img_it("./plots/xi_hat_%1.png"     ,f_xi_hat      ,res,512,512,-30,-30,30,30,mode);
    plot_img_it("./plots/xi_form0_%1.png"   ,f_xi_form0    ,res,512,512,-30,-30,30,30,mode);
    plot_img_it("./plots/xi_form1_%1.png"   ,f_xi_form1    ,res,512,512,-30,-30,30,30,mode);
    plot_img_it("./plots/gamma_%1.png"      ,f_gamma       ,res,512,512,-10,-10,10,10,mode);
    plot_img_it("./plots/digamma_%1.png"    ,f_digamma     ,res,512,512,-10,-10,10,10,mode);
    plot_img_it("./plots/beta_%1.png"       ,f_rbeta       ,res,512,512,-5,-5,5,5,mode);

    plot_img_it("./plots/tmoebius_zeta_form0_%1.png",[](cplx z,int it)->cplx{return f_zeta_form0(t_mobius(z,-1,1,1,1),it);},res,512,512,  -1.5,-0.5,-0.5,0.5,mode);
    plot_img_it("./plots/tmoebius_zeta_form1_%1.png",[](cplx z,int it)->cplx{return f_zeta_form1(t_mobius(z,-1,1,1,1),it);},res,512,512,  -1.5,-0.5,-0.5,0.5,mode);
    plot_img_it("./plots/tmoebius_zeta_form2_%1.png",[](cplx z,int it)->cplx{return f_zeta_form2(t_mobius(z,-1,1,1,1),it);},res,512,512,  -1.5,-0.5,-0.5,0.5,mode);
    plot_img_it("./plots/tmoebius_zeta_form3_%1.png",[](cplx z,int it)->cplx{return f_zeta_form3(t_mobius(z,-1,1,1,1),it);},res,512,512,  -1.5,-0.5,-0.5,0.5,mode);

    //DL ---------------------------------------------------------------------
    plot_img_it("./plots/TL_poles_1_2_%1.png"   ,TL_poles_1_2   ,res,64,64,-3,-3,3,3,mode);
    plot_img_it("./plots/T_poles_1_2_%1.png"    ,T_poles_1_2    ,res,64,64,-3,-3,3,3,mode);
    plot_img_it("./plots/L_poles_1_2_%1.png"    ,L_poles_1_2    ,res,64,64,-3,-3,3,3,mode);
    plot_img_it("./plots/geomTL_%1.png"         ,TorL_geometric ,res,256,256,-2,-2,2,2,mode);
    plot_img_it("./plots/geomL_%1.png"          ,L_geometric    ,res,2,128,-2,-2,2,2,mode);
    plot_img_it("./plots/geomT_%1.png"          ,T_geometric    ,res,2,64,-2,-2,2,2,mode);
    plot_img_it("./plots/expT_%1.png"           ,T_exp          ,res,2,64,-10,-10,10,10,mode);
    plot_img_it("./plots/sinT_%1.png"           ,T_sin          ,res,2,64,-10,-10,10,10,mode);
    plot_img_it("./plots/sinhT_%1.png"          ,T_sinh         ,res,2,64,-10,-10,10,10,mode);
    plot_img_it("./plots/cosT_%1.png"           ,T_cos          ,res,2,64,-10,-10,10,10,mode);
    plot_img_it("./plots/coshT_%1.png"          ,T_cosh         ,res,2,64,-10,-10,10,10,mode);
    plot_img_it("./plots/tanT_%1.png"           ,T_tan          ,res,2,64,-3,-3,3,3,mode);
    plot_img_it("./plots/tanhT_%1.png"          ,T_tanh         ,res,2,64,-3,-3,3,3,mode);
    plot_img_it("./plots/cotT_%1.png"           ,T_cot          ,res,2,64,-4,-4,4,4,mode);
    plot_img_it("./plots/cothT_%1.png"          ,T_coth         ,res,2,64,-4,-4,4,4,mode);
    plot_img_it("./plots/lnT_%1.png"            ,T_ln           ,res,2,64,-2,-2,2,2,mode);
    plot_img_it("./plots/bernoulliT_%1.png"     ,T_bernoulli    ,res,2,64,-10,-10,10,10,mode);
    plot_img_it("./plots/sqrtT_%1.png"          ,T_sqrt         ,res,2,64,-2,-2,2,2,mode);



    //MISC
    plot_img_it("./plots/TT_geomertic_ac_%1.png",TT_geomertic_ac,res,40,40,-2,-2,2,2,mode);
    plot_img_it("./plots/f_wallis_%1.png" , f_wallis ,res,512,512,-20,-20,20,20,MODE_ARG_ONLY);
}

void plot_nt()
{
    RealFuncList list1;
    list1.push_back( f_pi );
    list1.push_back( f_li_d );
    list1.push_back( f_r_d<100> );
    plot_real("./plots/f_pi100.csv","x;f_isprime;f_mu;f_pi;f_li;f_r<100>",list1,1000,1,100);

    RealFuncList list2;
    list2.push_back( f_pi );
    list2.push_back( f_r_d<100,0> );
    list2.push_back( f_r_d<100,2> );
    list2.push_back( f_r_d<100,4> );
    list2.push_back( f_r_d<100,8> );
    list2.push_back( f_r_d<100,16> );
    list2.push_back( f_r_d<100,32> );
    list2.push_back( f_r_d<100,512> );
    plot_real("./plots/f_r100.csv","x;f_pi;f_r<100,0>;f_r<100,2>;f_r<100,4>;f_r<100,8>;f_r<100,16>;f_r<350,32>",list2,1000,2,50);
}

cplx solve(std::function<cplx (cplx)> f,std::function<cplx (cplx)> fd,cplx w0,cplx z)
{
    //f(w)=z
    cplx w=w0;
    for(int i=0;i<20;i++){w=w+(z-f(w))/fd(w);}
    return w;
}

void plotRiemannSurfaceSqrtN(double R,double dr,double dt,int n,const char * filename)
{
    std::cout<<filename<<std::endl;

    std::ofstream os(filename);
    do
    {
        double theta=0.0;
        cplx w0=cplx(std::pow(R,1.0/n),0),w=w0;
        do
        {
            theta+=dt;
            cplx z( R*cos(theta), R*sin(theta) );
            w=solve([=](cplx s)->cplx{ return std::pow(s,n); },[=](cplx s)->cplx{ return double(n)*std::pow(s,n-1.0); },w,z);

            os<<z.real()<<";"<<z.imag()<<";"<<w.real()<<";" <<w.imag()<<";"<<std::abs(w)<<";"<<std::arg(w)<<std::endl;
        }
        while(std::abs(w-w0)>1e-5);
        R-=dr;
    }
    while(R>0);

    os.close();
}

void plotRiemannSurfaceLog(double R,double dr,double dt,int n,const char * filename)
{
    std::cout<<filename<<std::endl;

    std::ofstream os(filename);
    do
    {
        double theta=0.0;
        cplx w0=cplx(log(R),0),w=w0;
        do
        {
            theta+=dt;
            cplx z( R*cos(theta), R*sin(theta) );
            w=solve([=](cplx s)->cplx{ return std::exp(s); },[=](cplx s)->cplx{ return std::exp(s); },w,z);

            os<<z.real()<<";"<<z.imag()<<";"<<w.real()<<";" <<w.imag()<<";"<<std::abs(w)<<";"<<std::arg(w)<<std::endl;
        }
        while(std::abs(w-w0)>1e-5 && theta<n*2*M_PI);
        R-=dr;
    }
    while(R>0);

    os.close();
}


int main(int argc, char *argv[])
{
    Q_UNUSED(argc);
    Q_UNUSED(argv);

    std::cout<<sizeof(TYPE)<<" (eps,min,max)=("<<TYPE_EPSILON<<","<<TYPE_MIN<<","<<TYPE_MAX<<")"<<std::endl;

//    plotRiemannSurfaceSqrtN(2,0.01,2*M_PI/1000,2,"./plots/riemann_surface_sqrt_2.csv");
//    plotRiemannSurfaceSqrtN(2,0.01,2*M_PI/1000,3,"./plots/riemann_surface_sqrt_3.csv");
//    plotRiemannSurfaceSqrtN(2,0.01,2*M_PI/1000,4,"./plots/riemann_surface_sqrt_4.csv");
//    plotRiemannSurfaceSqrtN(2,0.01,2*M_PI/1000,5,"./plots/riemann_surface_sqrt_5.csv");
//    plotRiemannSurfaceLog  (5,0.1,2*M_PI/100,3,"./plots/riemann_log.csv");

    //bench_numbers();
    //bench_functions();
    //bench_TL_functions();

    //plot_nt();

    plot_img_it("./plots/erfT_%1.png"          ,T_erf           ,Resolution(1024),2,128,-5,-5,5,5,MODE_ARG_ONLY);
    plot_img_it("./plots/CT_%1.png"            ,T_C             ,Resolution(1024),2,64,-5,-5,5,5,MODE_ARG_ONLY);
    plot_img_it("./plots/ST_%1.png"            ,T_S             ,Resolution(1024),2,64,-5,-5,5,5,MODE_ARG_ONLY);

    //plot_all(Resolution(1024),MODE_ARG_ONLY);
    //plot_fractals(Resolution(1024));

    //plot_img("./plots/sqrt.png"     ,func_sqrt      ,Resolution(512),-2,-2,2,2,MODE_ARG_AND_MOD_ISO);

    //plot_img_it("./plots/f_zeta_prime_form0_k1_%1.png" , f_zeta_prime_form0 ,Resolution(512,512),512,512,-20,-20,20,20,MODE_ARG_ONLY);
    //plot_img_it("./plots/f_zeta_prime_form1_k1_%1.png" , f_zeta_prime_form1 ,Resolution(512,512),512,512,-20,-20,20,20,MODE_ARG_ONLY);
    //plot_img_it("./plots/beta_%1.png"       ,f_rbeta       ,256,512,512,-5,-5,5,5,MODE_BUMP);
    //plot_img_it("./plots/zeta_form4_%1.png" ,f_zeta_form4  ,512,4096,4096,-30,-30,30,30,MODE_MOD_ONLY);
    //plot_img_it("./plots/zeta_form4_%1.png" ,f_zeta_form5  ,Resolution(200,1600),1024,1024,-10,-80,10,80,MODE_BUMP);
    //plot_img_it("./plots/zeta_form4_%1.png" ,f_zeta_form5  ,Resolution(200,1600),1024,1024,-10,-80,10,80,MODE_ARG_ONLY);

    //bench_hypercplx();

    return 0;
}
