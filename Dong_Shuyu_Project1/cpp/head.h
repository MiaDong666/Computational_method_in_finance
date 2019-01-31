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