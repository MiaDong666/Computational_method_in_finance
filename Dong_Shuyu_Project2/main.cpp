#include <iostream>
#include "head.h"
#include <math.h>
#include<algorithm>
#include<fstream>
#include<sstream>
#include<windows.h>
#include <ctime>
#include <stdlib.h>
#include<stdio.h>

using namespace std;

int main() {
    //Q1
    int n=1000;
    double vx,vy,co,meanx,meany;
    vx=1;
    vy=1;
    co=-0.7;
    meanx=0;
    meany=0;
    pair<double*,double*> pair1=Bi_Normal(n,vx,vy,co,meanx,meany);
    double*X=pair1.first;
    double*Y=pair1.second;
    double rho;
    rho=corr(n,X,Y);
    cout<<"Q1:"<<endl;
    cout<<"Rho is: "<<rho<<endl;
    cout<<endl;


    //Q2
    //generate random numbers
    n=10000;
    vx=1;
    vy=1;
    co=0.6;
    meanx=0;
    meany=0;
    pair<double*,double*> pair2=Bi_Normal(n,vx,vy,co,meanx,meany);
    X=pair2.first;
    Y=pair2.second;
    //monte-carlo simulation
    int N=10000;
    double *E2=new double[N];
    for (int i = 0; i < N; ++i)
    {
        double t=X[i]*X[i]*X[i]+sin(Y[i])+X[i]*X[i]*Y[i];
        if(t>0)
        {
            E2[i]=t;
        }
        else{
            E2[i]=0;
        }
    }
    cout<<"Q2:"<<endl;
    cout<<"Value of E is: "<<mean(E2,N)<<endl;
    cout<<endl;

    //Q3
    //a}Ea1
    N=10000;
    double *Ea1=new double[N];
    double*W5=Wiener(N,5);

    for(int i=0;i<N;i++)
    {
        Ea1[i]=W5[i]*W5[i]+sin(W5[i]);
    }
    cout<<"Q3 a):"<<endl;
    cout<<"Expectation of Ea1 is: "<<mean(Ea1,N)<<". Variance is: "<<var(Ea1,N)<<endl;
    //Ea2
    double t;
    t=0.5;
    double*Wt2=Wiener(N,t);
    double*Ea2=new double[N];
    for(int i=0;i<N;i++)
    {
        Ea2[i]=exp(t/2)*cos(Wt2[i]);
    }
    cout<<"Expectation of Ea2 is: "<<mean(Ea2,N)<<". Variance is: "<<var(Ea2,N)<<endl;
    //Ea3
    t=3.2;
    double*Wt3=Wiener(N,t);
    double*Ea3=new double[N];
    for(int i=0;i<N;i++)
    {
        Ea3[i]=exp(t/2)*cos(Wt3[i]);
    }
    cout<<"Expectation of Ea3 is: "<<mean(Ea3,N)<<". Variance is: "<<var(Ea3,N)<<endl;
    //Ea4
    t=6.5;
    double*Wt4=Wiener(N,t);
    double*Ea4=new double[N];
    for(int i=0;i<N;i++)
    {
        Ea4[i]=exp(t/2)*cos(Wt4[i]);
    }
    cout<<"Expectation of Ea4 is: "<<mean(Ea4,N)<<". Variance is: "<<var(Ea4,N)<<endl;



    //Q3 c)Eb1
    //take Y=sin(Wt), delta=0
    double gam;
    double* Y_Q3=new double[N];
    for (int i=0;i<N;i++)
    {
        Y_Q3[i]=W5[i]*W5[i];
    }
    cout<<"Q3 c):"<<endl;
    cout<<"By controlling variate using Y=Wt*Wt,"<<endl;
    double *Test=controlvar(Ea1,Y_Q3,5,N);
    cout<<"Expectation of Eb1 is: "<<mean(Test,N)<<". Variance is: "<<var(Test,N)<<endl;

    //Q3 c)Ea2
    for (int i=0;i<N;i++)
    {
        Y_Q3[i]=Wt2[i]*Wt2[i];
    }
    Test=controlvar(Ea2,Y_Q3,0.5,N);
    cout<<"Expectation of Eb2 is: "<<mean(Test,N)<<
                                   ". Variance is: "<<var(Test,N)<<endl;
    //Q3 c)Ea3
    for (int i=0;i<N;i++)
    {
        Y_Q3[i]=Wt3[i]*Wt3[i];
    }
    Test=controlvar(Ea3,Y_Q3,3.2,N);
    cout<<"Expectation of Eb3 is: "<<mean(Test,N)<<
        ". Variance is: "<<var(Test,N)<<endl;
    //Q3 c)Ea4
    for (int i=0;i<N;i++)
    {
        Y_Q3[i]=Wt4[i]*Wt4[i];
    }
    Test=controlvar(Ea4,Y_Q3,6.5,N);
    cout<<"Expectation of Eb4 is: "<<mean(Test,N)<<
                                   ". Variance is: "<<var(Test,N)<<endl;
    cout<<endl;




    //Q4 a)
    t=5;
    N=10000;
    double*Wt=Wiener(N,t);
    double sigma=0.2;
    double r=0.04;
    double S0=88;
    double*Sto1=Stock(sigma,r,t,Wt,S0,N);
    double *call1=new double[N];
    double K=100;
    for (int i = 0; i < N; i++) {
        call1[i]=exp(-r*t)*Callpayoff(Sto1[i],K);
    }
    cout<<"Q4："<<endl;
    cout<<"The price of call option is: "<<mean(call1,N)<<endl;

    //b)
    cout<<"The Black-Sholes formula option price is: "<<Black_Sholes(sigma,r,t,S0,K)<<endl;
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
    cout<<"Using variance-reduction method, we get the price of the call option is: "<<mean(call_v,N)<<endl;
    cout<<endl;

    //Q5 a)
    N=1000;
    r=0.04;
    sigma=0.18;
    S0=88;
    double*W_1=Wiener(N,1);
    double*S_1=Stock(sigma,r,1,W_1,S0,N);
    double*ES=new double[10];
    ES[0]=mean(S_1,N);
    double*W_2=Wiener(N,2);
    double*S_2=Stock(sigma,r,2,W_2,S0,N);
    ES[1]=mean(S_2,N);
    double*W_3=Wiener(N,3);
    double*S_3=Stock(sigma,r,3,W_3,S0,N);
    ES[2]=mean(S_3,N);
    double*W_4=Wiener(N,4);
    double*S_4=Stock(sigma,r,4,W_4,S0,N);
    ES[3]=mean(S_4,N);
    double*W_5=Wiener(N,5);
    double*S_5=Stock(sigma,r,5,W_5,S0,N);
    ES[4]=mean(S_5,N);
    double*W_6=Wiener(N,6);
    double*S_6=Stock(sigma,r,6,W_6,S0,N);
    ES[5]=mean(S_6,N);
    double*W_7=Wiener(N,7);
    double*S_7=Stock(sigma,r,7,W_7,S0,N);
    ES[6]=mean(S_7,N);
    double*W_8=Wiener(N,8);
    double*S_8=Stock(sigma,r,8,W_8,S0,N);
    ES[7]=mean(S_8,N);
    double*W_9=Wiener(N,9);
    double*S_9=Stock(sigma,r,9,W_9,S0,N);
    ES[8]=mean(S_9,N);
    double*W_10=Wiener(N,10);
    double*S_10=Stock(sigma,r,10,W_10,S0,N);
    ES[9]=mean(S_10,N);


    //output for graph
    ofstream file1;
    file1.open("C:\\Users\\Mia Dong\\Desktop\\Dong_Shuyu_Project2\\Shuyu_Stock1.csv");
    for (int i = 0; i < 10; i++)
    {
        file1 << ES[i] <<  ",\n";
    }
    file1 << "\n";
    file1.close();

    //Q5 b)
    N=1000;
    t=10;
    double *path_1=stocksimulation(sigma,r,t,S0,N);
    double *path_2=stocksimulation(sigma,r,t,S0,N);
    double *path_3=stocksimulation(sigma,r,t,S0,N);
    double *path_4=stocksimulation(sigma,r,t,S0,N);
    double *path_5=stocksimulation(sigma,r,t,S0,N);
    double *path_6=stocksimulation(sigma,r,t,S0,N);

    //output for graph
    ofstream file2;
    file2.open("C:\\Users\\Mia Dong\\Desktop\\Dong_Shuyu_Project2\\Shuyu_Stockpath1.csv");
    for (int i = 0; i < N+1; i++)
    {
        file2 << path_1[i]<<","<<path_2[i]<<","<<path_3[i]<<","<<path_4[i]<<","<<path_5[i]<<","<<path_6[i] <<  ",\n";
    }
    file2 << "\n";
    file2.close();

    //Q5 c)
    sigma=0.35;
    S_1=Stock(sigma,r,1,W_1,S0,N);
    double*ES_2=new double[10];
    ES_2[0]=mean(S_1,N);
    S_2=Stock(sigma,r,2,W_2,S0,N);
    ES_2[1]=mean(S_2,N);
    S_3=Stock(sigma,r,3,W_3,S0,N);
    ES_2[2]=mean(S_3,N);
    S_4=Stock(sigma,r,4,W_4,S0,N);
    ES_2[3]=mean(S_4,N);
    S_5=Stock(sigma,r,5,W_5,S0,N);
    ES_2[4]=mean(S_5,N);
    S_6=Stock(sigma,r,6,W_6,S0,N);
    ES_2[5]=mean(S_6,N);
    S_7=Stock(sigma,r,7,W_7,S0,N);
    ES_2[6]=mean(S_7,N);
    S_8=Stock(sigma,r,8,W_8,S0,N);
    ES_2[7]=mean(S_8,N);
    S_9=Stock(sigma,r,9,W_9,S0,N);
    ES_2[8]=mean(S_9,N);
    S_10=Stock(sigma,r,10,W_10,S0,N);
    ES_2[9]=mean(S_10,N);


    //output for graph
    ofstream file3;
    file3.open("C:\\Users\\Mia Dong\\Desktop\\Dong_Shuyu_Project2\\Shuyu_Stock2.csv");
    for (int i = 0; i < 10; i++)
    {
        file3 << ES[i] <<  ",\n";
    }
    file3 << "\n";
    file3.close();

    //Q5 c)2
    path_1=stocksimulation(sigma,r,t,S0,N);
    path_2=stocksimulation(sigma,r,t,S0,N);
    path_3=stocksimulation(sigma,r,t,S0,N);
    path_4=stocksimulation(sigma,r,t,S0,N);
    path_5=stocksimulation(sigma,r,t,S0,N);
    path_6=stocksimulation(sigma,r,t,S0,N);

    //output for graph
    ofstream file4;
    file4.open("C:\\Users\\Mia Dong\\Desktop\\Dong_Shuyu_Project2\\Shuyu_Stockpath2.csv");
    for (int i = 0; i < N+1; i++)
    {
        file4 << path_1[i]<<","<<path_2[i]<<","<<path_3[i]<<","<<path_4[i]<<","<<path_5[i]<<","<<path_6[i] <<  ",\n";
    }
    file4 << "\n";
    file4.close();

    //Q6 a)Euler scheme
    N=10000;
    double interval=1;
    double dx;
    dx=interval/N;
    double sum_Euler;
    sum_Euler=0;
    double fx;
    double*X0=new double[N+1];
    X0[0]=0;
    for (int i = 0; i < N; ++i) {

        fx=sqrt(1-X0[i]*X0[i]);
        sum_Euler=sum_Euler+dx*fx;
        X0[i+1]=X0[i]+dx;
    }
    cout<<"Q6:"<<endl;
    cout<<"Euler's scheme gives that pi is: "<<4*sum_Euler<<endl;

    //Q6 b)Monte-Carlo
    double* fu=new double[N];
    double* Uni=get_uniform(N);
    for (int i = 0; i <N ; i++) {
        fu[i]=sqrt(1-Uni[i]*Uni[i]);
    }

    cout<<"Monte-Carlo simulation gives that pi is: "<<4*mean(fu,N)<<". Variance is: "<<var(fu,N)<<endl;


    //Q6 c)
    //using acceptance-rejection to generate Y;take g(x)=1.5
    double*Y_test=get_uniform(2*N);
    double*X_test=new double[N];
    double*U_test=get_uniform(2*N);
    int count=0;
    double alpha=0.74;
    double *fy=new double[2*N];
    double *gy=new double[2*N];
    double *hy=new double[2*N];
    for(int i=0;i<2*N;i++){
        fy[i]=(1-alpha*Y_test[i]*Y_test[i])/(1-alpha/3);
        gy[i]=1.5;
        hy[i]=fy[i]/gy[i];
        if(U_test[i]<hy[i]){
            X_test[count]=Y_test[i];
            count=count+1;
        }
        if(count>=N){
            break;
        }
    }

    double* hu=new double[N];
    double* gu=new double[N];

    for (int i=0;i<N;i++){
        hu[i]=(1-alpha*X_test[i]*X_test[i])/(1-alpha/3);
        fu[i]=sqrt(1-X_test[i]*X_test[i]);
        gu[i]=fu[i]/hu[i];
    }
    cout<<"Using importance sampling, pi is: "<<4*mean(gu,N)<<". Variance is: "<<var(gu,N)<<endl;








    return 0;
}