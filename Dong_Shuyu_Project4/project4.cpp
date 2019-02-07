//
// Created by Mia Dong on 2/3/2019.
//
#include "head.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

//type=0,a);type=1,b);type=2,c);type=3,d)
double Binomial_Europeanoption(int N,double sigma,double r, double t,double s0,double K,int type){
    double delta=t/N;
    double d;
    double u;
    double p;
    switch (type){
        case 0:{
            double c=0.5*(exp(-r*delta)+exp((r+sigma*sigma)*delta));
            d=c-sqrt(c*c-1);
            u=1.0/d;
            p=(exp(r*delta)-d)/(u-d);
            break;
        }
        case 1:{
            p=0.5;
            d=exp(r*delta)*(1-sqrt(exp(sigma*sigma* delta)-1));
            u=exp(r*delta)*(1+sqrt(exp(sigma*sigma* delta)-1));
            break;
        }
        case 2:{
            p=0.5;
            u=exp((r-0.5*sigma*sigma)* delta+sigma*sqrt(delta));
            d=exp((r-0.5*sigma*sigma)* delta-sigma*sqrt(delta));
            break;
        }
        default:{
            u=exp(sigma*sqrt(delta));
            d=exp(-sigma*sqrt(delta));
            p=0.5+0.5*(r-0.5*sigma*sigma)*sqrt(delta)/sigma;
            break;
        }
    }
   // cout<<p<<endl;
    //creating nodes
    double uu=u/d;
    double R=exp(r*delta);
    vector<double>price(N+1);
    price[0]=s0*pow(d,N);

    for(int i=1;i<(N+1);i++){
        price[i]=price[i-1]*uu;
    }

    vector<double>call_value(N+1);
    for(int i=0;i<(N+1);i++){
        call_value[i]=Callpayoff(price[i],K);
    }
    for(N=N-1;N>=0;--N){
        for(int i=0;i<(N+1);i++){
            call_value[i]=(p*call_value[i+1]+(1-p)*call_value[i])/R;
        }
    }
    return call_value[0];
};

//type=0, European option; type=1, American
double Binomial_put(int N,double sigma,double r, double t,double s0,double K,int type){
    double delta=t/N;
    double u=exp(sigma*sqrt(delta));
    double d=exp(-sigma*sqrt(delta));
    double p=0.5+0.5*(r-0.5*sigma*sigma)*sqrt(delta)/sigma;

    //creating nodes
    double dd=d/u;
    double R=exp(r*delta);
    double**price=new double*[N+1];
    double** put_value=new double*[N+1];
    for(int i=N;i>=0;--i){
        price[i]=new double[N+1];
        put_value[i]=new double[N+1];
    }


    for(int i=0;i<(N+1);i++){
        price[i][0]=s0*pow(u,i);
        for(int j=1;j<i+1;j++){
            price[i][j]=price[i][j-1]*dd;
        }
    }


    for(int i=0;i<(N+1);i++){
        put_value[N][i]=max(0.0,K-price[N][i]);
    }
    for(int i=N-1;i>=0;--i){
        for(int j=0;j<(i+1);j++){
            switch(type){
                case 0:{
                    put_value[i][j]=(p*put_value[i+1][j]+(1-p)*put_value[i+1][j+1])/R;
                    break;
                }
                default:{
                    put_value[i][j]=max((K-price[i][j]),((p*put_value[i+1][j]+(1-p)*put_value[i+1][j+1])/R));
                    break;
                }
            }
        }
    }

    return put_value[0][0];
}

