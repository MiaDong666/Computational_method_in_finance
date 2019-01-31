//
// Created by Mia Dong on 1/21/2019.
//

#ifndef DONG_SHUYU_PROJECT2_HEAD_H
#define DONG_SHUYU_PROJECT2_HEAD_H
#include<iostream>
#include <utility>
#include <string>
using namespace std;
#endif //DONG_SHUYU_PROJECT2_HEAD_H
#pragma once
double* get_uniform(int n);
double* get_uniform_UN(int n);
double mean(double* X,int n);
double var(double* X, int n);
double* bernoulli(double p,int n);
double* Binom(double p, int k, int n);
double* Exp(double lambda,int n);
double probbig(double * X,int n, double rule);
double* normal_BM(int n);
double* normal_PM(int n);
pair<double*,double*> Bi_Normal(int n, double vx,double vy,double cov,double meanx,double meany);
double cov(int n, double *X,double *Y);
double corr(int n, double*X,double*Y);
double* Wiener(int n, double T);
double* Stock(double sigma,double r,double t,double* W,double s0,int N);
double Callpayoff(double S,double K);
double* controlvar(double*X,double*Y,double delta,int N);
double* stocksimulation(double sigma,double r,double time,double S0, int N);
double pnorm(double x);
double Black_Sholes(double sigma,double r,double t,double S0,double K);