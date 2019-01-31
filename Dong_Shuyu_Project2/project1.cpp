//
// Created by Mia Dong on 1/21/2019.
//

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
#include<cstdio>

using namespace std;
using namespace std::chrono;

//LGM-generator
double* get_uniform(int n)
{
    double* UN= new double[n];
    long long * U = new long long[n];
    long long a, b, m;
    static long x0=clock();
    a = pow(7,5);
    b = 0;
    m = pow(2,31)-1;
    x0=(x0*231+937)%10007;
    U[0] = x0;
    for (int i=1;i<=n;i++)
    {
        U[i] = (U[i - 1] * a + b) %m;
        int x = 0;
    }
    for (int i = 0; i < n; i++)
    {
        double c  = ((double)U[i])/((double)m-1);
        UN[i] = c;
    }


    return UN;
}

double* get_uniform_UN(int n)
{
    double* UN = new double[n];
    for (int i = 0; i < n; i++)
    {
        UN[i] = rand()*1.0/(RAND_MAX*1.0);
    }

    return UN;
}
double mean(double * X,int n)
{
    double sum, mean;
    sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum = sum + X[i];
    }
    mean = sum / n;
    return mean;
}
double var(double* X, int n)
{
    double sum,variance;
    sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum = sum + (X[i] - mean(X, n))*(X[i] - mean(X, n));
    }
    variance = sum / (n - 1);
    return variance;
}
double* bernoulli(double p,int n)
{
    double* UN;
    double* X = new double[n];
    UN = get_uniform(n);
    for (int i = 0; i < n; i++)
    {
        if (UN[i] <= p)
        {
            X[i] = 1;
        }
        else
        {
            X[i]= 0;
        }

    }
    return X;
}

double* Binom(double p, int k, int n)
{
    int N = k * n;
    double * Ber = bernoulli(p,N);
    double * Bi = new double[n];
    for (int i=0;i<n;i++)
    {
        double sum;
        sum = 0;
        for (int j=k*i;j<k*(i+1);j++)
        {
            sum = sum + Ber[j];
        }
        Bi[i] = sum;
    }


    return Bi;
}
double* Exp(double lambda,int n)
{
    double* U = get_uniform(n);
    double* X = new double[n];


    for (int i=0;i<n;i++)
    {
        X[i] = -lambda * log(U[i]);
    }
    return X;

}
//probability if X>= certain number
double probbig(double * X,int n, double rule)
{
    double p;
    int count = 0;
    for (int i=0;i<n;i++)
    {
        if (X[i]>=rule)
        {
            count = count + 1;
        }
    }
    p = (double)count / (double)n;
    return p;
}
//Box-Muller method
double* normal_BM(int n)
{
    double *Normal=new double[n];
    double *U;
    U = get_uniform(n);
    for (int i=0;i<(n/2);i++)
    {
        Normal[2 * i] = sqrt(-2 * log(U[2 * i]))*cos(2 * M_PI*U[2 * i + 1]);
        Normal[2*i+1]= sqrt(-2 * log(U[2 * i]))*sin(2 * M_PI*U[2 * i +1 ]);
    }
    return Normal;
}
//Polar-Masarglia method
double* normal_PM(int n)
{
    double *Normal=new double[n];
    double *U=get_uniform(2*n);
    int count=0;
    double V1,V2,W;
    for(int i=0;i <n;i++)
    {
        V1=U[2*i]*2-1;
        V2=U[2*i+1]*2-1;
        W=V1*V1+V2*V2;
        if(W<=1)
        {

            Normal[2*count]=V1*sqrt(-2*log(W)/W);
            Normal[2*count+1]=V2*sqrt(-2*log(W)/W);
            count=count+1;
        }
        if (count>=n/2)
        {
            break;
        }
    }
    return Normal;
}