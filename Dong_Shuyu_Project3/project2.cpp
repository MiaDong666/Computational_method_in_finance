//
// Created by Mia Dong on 1/21/2019.
//

// random number generator
#include<ctime>
#include<iostream>
#include<math.h>
#include<random>
#include<stdio.h>
#include"head.h"
#include <algorithm>
#include <chrono>
#include <utility>
#include <string>
#include <cmath>

using namespace std;
using namespace std::chrono;

//Bivariate-Normal
pair<double*,double*> Bi_Normal(int n, double vx,double vy,double cov, double meanx,double meany)
{
    double rho,sigmax,sigmay;
    double*X=new double[n];
    double*Y=new double[n];
    double*Z1=normal_BM(n);
    double*Z2=normal_PM(n);
    rho=cov/sqrt(vx*vy);
    sigmax=sqrt(vx);
    sigmay=sqrt(vy);
    for(int i=0;i<n;i++)
    {
        X[i]=meanx+sigmax*Z1[i];
        Y[i]=meany+sigmay*rho*Z1[i]+sqrt(1-rho*rho)*sigmay*Z2[i];
    }
    pair<double*,double*> result(X,Y);
    return result;

}
//cov function
double cov(int n, double *X,double *Y)
{
    double sum,cov;
    sum=0;
    for(int i=0;i<n;i++)
    {
        sum=sum+(X[i]-mean(X,n))*(Y[i]-mean(Y,n));
    }
    cov=sum/(n-1);
    return cov;
}

//correlation function
double corr(int n, double*X,double*Y)
{
    double rho;
    rho=cov(n,X,Y)/sqrt(var(X,n)*var(Y,n));
    return rho;
};
//Wiener process;
double* Wiener(int n,double T)
{
    double*Normal=normal_BM(n);
    double *W = new double[n];
    for(int i=0;i<n;i++)
    {
        W[i]=sqrt(T)*Normal[i];
    }
    return W;
};
double* Stock(double sigma,double r,double t,double* W,double S0,int N){
    double* S=new double[N];
    for(int i=0;i<N;i++){
        S[i]=S0*exp(sigma*W[i]+ (r-0.5*sigma*sigma)*t);
    }
    return S;
}
double Callpayoff(double S,double K)
{
    double c;
    if(S>K)
    {
        c=S-K;
    } else{
        c=0;
    }
    return c;
}
double *controlvar(double*X,double*Y,double delta,int N){
    double gam,meanT;
    double*T= new double[N];
    gam=cov(N,X,Y)/var(Y,N);
    for (int i = 0; i <N ; i++) {
        T[i]=X[i]-gam*(Y[i]-delta);
    }
    return T;
}
double*stocksimulation(double sigma,double r,double time,double S0, int N ){
    double* path=new double[N+1];
    double dt=time/N;
    double*Normal=normal_BM(N);
    path[0]=S0;
    double ds;
    for (int i = 1; i < N+1; ++i) {
        ds=r*path[i-1]*dt+sigma*path[i-1]*sqrt(dt)*Normal[i-1];
        path[i]=path[i-1]+ds;
    }

    return path;
}


double Black_Sholes(double sigma,double r,double t,double S0,double K){
    double d1=1/(sigma*sqrt(t))*(log(S0/K)+(r+0.5*sigma*sigma)*t);
    double d2=d1-sigma*sqrt(t);
    double c=pnorm(d1)*S0-pnorm(d2)*K*exp(-r*t);
    return c;
}




