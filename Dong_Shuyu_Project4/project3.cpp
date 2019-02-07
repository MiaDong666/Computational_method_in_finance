//
// Created by Mia Dong on 1/28/2019.
//
#include "head.h"
#include <math.h>
#include <algorithm>
#include <utility>
#include <chrono>
#include <string>
using namespace std;
using namespace std::chrono;

double callopition_montecarlo(double sigma,double r,double t,double S0,double K){
    int N=10000;
    double*Wt=Wiener(N,t);
    double*Sto1=Stock(sigma,r,t,Wt,S0,N);
    double *call1=new double[N];
    for (int i = 0; i < N; i++) {
        call1[i]=exp(-r*t)*Callpayoff(Sto1[i],K);
    }
    double*Wt_=new double[N];
    for (int i = 0; i < N; i++) {
        Wt_[i]=-Wt[i];
    }
    double*Sto2=Stock(sigma,r,t,Wt_,S0,N);
    double*call2=new double[N];
    for (int i = 0; i < N; i++) {
        call2[i]=exp(-r*t)*Callpayoff(Sto2[i],K);
    }
    double*call_v=new double[N];
    for (int i = 0; i < N; i++) {
        call_v[i]=0.5*(call1[i]+call2[i]);
    }

    return mean(call_v,N);
};
//normal distribution cdf
long double Nx(double x){
    long double Nx;
    long double d1=0.0498673470;
    long double d2 = 0.0211410061;
    long double d3 = 0.0032776263;
    long double d4 = 0.0000380036;
    long double d5 = 0.0000488906;
    long double d6 = 0.0000053830;
    long double temp=1+d1*x+d2*x*x+d3*x*x*x+d4*pow(x,4)+d5*pow(x,5)+d6*pow(x,6);
    Nx=1-0.5*pow(temp,-16);
    return Nx;
}
double pnorm(double x)
{

    double norm;
    if(x>=0){
        norm=Nx(x);
    }else{
        norm=1-Nx(-x);
    }
    return norm;
}
double delta(double sigma,double r,double t,double S0,double K){
    double episilon=0.5;
    double Su=callopition_montecarlo(sigma,r,t,(S0+episilon),K);
    double Sd=callopition_montecarlo(sigma,r,t,(S0-episilon),K);
    double delta=(Su-Sd)/(2*episilon);
    return delta;
}
double theta(double sigma,double r,double t,double S0,double K){
    double episilon=0.3;
    double Su=callopition_montecarlo(sigma,r,(t+episilon),S0,K);
    double Sd=callopition_montecarlo(sigma,r,t-episilon,S0,K);
    double theta=-(Su-Sd)/(2*episilon);
    return theta;
}
double vega(double sigma,double r,double t,double S0,double K){
    double episilon=0.1;
    double Su=callopition_montecarlo(sigma+episilon,r,t,S0,K);
    double Sd=callopition_montecarlo(sigma-episilon,r,t,S0,K);
    double Vega=(Su-Sd)/(2*episilon*100);
    return Vega;
}
double gamma(double sigma,double r,double t,double S0,double K){
    double episilon=2;
    double Su=callopition_montecarlo(sigma,r,t,S0+episilon,K);
    double Sd=callopition_montecarlo(sigma,r,t,S0-episilon,K);
    double S=callopition_montecarlo(sigma,r,t,S0,K);
    double gamma=(Su-2*S+Sd)/(episilon*episilon);
    return gamma;
}
double rho(double sigma,double r,double t,double S0,double K){
    double episilon=0.05;
    double Su=callopition_montecarlo(sigma,r+episilon,t,S0,K);
    double Sd=callopition_montecarlo(sigma,r-episilon,t,S0,K);
    double rho=(Su-Sd)/(2*episilon*100);
    return rho;
}

//Heston model,type=0-Full Truncation, type=1-Partial Truncation, type=2-Reflection

double Heston(double sigma,double r,double t,double S0,double V0,double alpha,double beta,double rho,int type){
    int N=1000;
    double*Vol=new double[N+1];
    double*S=new double[N+1];
    double dt=t/N;
    Vol[0]=V0;
    S[0]=S0;
    double dv,ds;
    //W1,W2
    pair<double*,double*> Wie=Bi_Normal(N,dt,dt,(rho*dt),0,0);
    double *W1=Wie.first;
    double *W2=Wie.second;
    //simulate Vt,St
    for (int i = 1; i < N+1; ++i) {
        if(Vol[i-1]<0){
            switch(type){
                //Full Truncation
                case 0:{
                    ds=r*S[i-1]*dt+sqrt(0)*S[i-1]*W1[i-1];
                    dv=alpha*(beta-0)*dt+sigma*sqrt(0)*W2[i-1];
                    Vol[i]=Vol[i-1]+dv;
                    S[i]=S[i-1]+ds;
                    break;
                }
                //Partial Truncation
                case 1:{
                    ds=r*S[i-1]*dt+sqrt(0)*S[i-1]*W1[i-1];
                    dv=alpha*(beta-Vol[i-1])*dt+sigma*sqrt(0)*W2[i-1];
                    Vol[i]=Vol[i-1]+dv;
                    S[i]=S[i-1]+ds;
                    break;
                }
                //reflection
                default:{
                    ds=r*S[i-1]*dt+sqrt(-Vol[i-1])*S[i-1]*W1[i-1];
                    dv=alpha*(beta+Vol[i-1])*dt+sigma*sqrt(-Vol[i-1])*W2[i-1];
                    Vol[i]=-Vol[i-1]+dv;
                    S[i]=S[i-1]+ds;
                    break;
                }
            }

        }else{
            ds=r*S[i-1]*dt+sqrt(Vol[i-1])*S[i-1]*W1[i-1];
            dv=alpha*(beta-Vol[i-1])*dt+sigma*sqrt(Vol[i-1])*W2[i-1];
            Vol[i]=Vol[i-1]+dv;
            S[i]=S[i-1]+ds;
        }
    }
    return S[N];
}
double* getHalton(int n,int base){
    double *r=new double[n];
    for(int i=0;i<n;i++){
        int k=i+1;
        double f=1;
        r[i]=0;
        do{
            f=f/base;
            r[i]=r[i]+f*(k%base);
            k=k/base;
        }while (k>0);
    }
    return r;
}
