#include <iostream>
#include "head.h"
#include <math.h>
#include <cmath>
#include<fstream>
#include<sstream>
using namespace std;
int main(){
    //Q1
    //M：times for monte-Carlo simulation
    int N=1000;
    int M=1000;

    double integral=2;
    double dt=integral/N;
    double* Mx=new double[M];

    //monte-carlo Xt
    for(int i=0;i<M;i++){
        //Xt
        double* W1=Wiener(N,dt);
        double*Xt=new double[N+1];
        Xt[0]= 1;
        double dx;
        for(int j=1;j<N+1;j++){
           dx=(0.2-0.5*Xt[j-1])*dt+2.0/3.0*W1[j-1];
           Xt[j]=Xt[j-1]+dx;
        }
        Mx[i]=Xt[N];
    }
    double* E1=new double[M];
    for(int i=0;i<M;i++){
        E1[i]=cbrt(Mx[i]);

    }
    cout<<"Q1："<<endl;
    cout<<"E1 (X^1/3) is: "<<mean(E1,M)<<endl;

    //monte-carlo Yt
    double* My=new double[M];
    for(int i=0;i<M;i++){
        //Yt
        double*W2=Wiener(N,dt);
        double*Yt=new double[N+1];
        Yt[0]=0.75;
        double t=0;
        double dy;
        for(int j=1;j<N+1;j++){
            dy=((2/(1+t))*Yt[j-1]+(1+pow(t,3))/3)*dt+(1+pow(t,3))/3*W2[j-1];
            t=t+dt;
            Yt[j]=Yt[j-1]+dy;
        }
        My[i]=Yt[N];
    }
    //probability
    cout<<"Pr(Y2>5) is: "<<probbig(My,M,5)<<endl;

    double* E3=new double[M];
    for(int i=0;i<M;i++){
        if(Mx[i]>1){
            E3[i]=Mx[i]*My[i];
        } else{
            E3[i]=0;
        }
    }
    cout<<"E3 (X2Y2) is: "<<mean(E3,M)<<endl;

    integral=3;
    dt=integral/N;
    for(int i=0;i<M;i++){
        //Yt
        double*W2=Wiener(N,dt);
        double*Yt=new double[N+1];
        Yt[0]=0.75;
        double t=0;
        double dy;
        for(int j=1;j<N+1;j++){
            dy=((2/(1+t))*Yt[j-1]+(1+pow(t,3))/3)*dt+(1+pow(t,3))/3*W2[j-1];
            t=t+dt;
            Yt[j]=Yt[j-1]+dy;
        }
        My[i]=Yt[N];
    }
    cout<<"E2 (Y3) is: "<<mean(My,M)<<endl;
    cout<<endl;

    //Q2
    //E1
    N=1000;
    M=1000;
    integral=3;
    dt=integral/N;
    double *E21=new double[M];

    for (int i = 0; i <M ; ++i) {
        double* W1=Wiener(N,dt);
        double* W2=Wiener(N,dt);
        double*Xt=new double[N+1];
        Xt[0]= 1;
        double dx;
        for(int j=1;j<N+1;j++){
            dx=0.25*Xt[j-1]*dt+1.0/3.0*Xt[j-1]*W1[j-1]-0.75*Xt[j-1]*W2[j-1];
            Xt[j]=Xt[j-1]+dx;
        }
        Mx[i]=Xt[N];
    }
    for(int i=0;i<M;i++){
        E21[i]=cbrt(1+Mx[i]);
    }
    cout<<"Q2: "<<endl;
    cout<<"E1 is: "<<mean(E21,M)<<endl;

    //E2
    integral=3;
    double* E22=new double[M];

    for (int i = 0; i <M ; ++i) {
        double* W1=Wiener(N,integral);
        double* W2=Wiener(N,integral);
        double*Yt=new double[N];
        Yt[i]=exp(-0.08*integral+1.0/3.0*W1[i]+0.75*W2[i]);
        My[i]=Yt[i];
    }
    for(int i=0;i<M;i++){
        E22[i]=cbrt(1+My[i]);
    }
    cout<<"E2 is: "<<mean(E22,M)<<endl;
    cout<<endl;


    //Q3
    //a&b
    double X=20,sigma=0.25,r=0.04,T=0.5;
    double S=20;
    cout<<"Q3: "<<endl;
    cout << "Given X=20, sigma=0.25, r=0.04, T=0.5, s0=20"<< endl;
    cout<<"Monte Carlo gives the call option price c1: "<<callopition_montecarlo(sigma,r,T,S,X)<<endl;
    cout<<"Black-Sholes gives the call option price c2: "<<Black_Sholes(sigma,r,T,S,X)<<endl;

    //c
    double* s0=new double[11];
    s0[0]=15;
    for(int i=1;i<11;i++){
        s0[i]=s0[i-1]+1;
    }
    double *del=new double[11];
    double *the=new double[11];
    double* veg=new double[11];
    double*gam=new double[11];
    double*rh=new double[11];
    double*cal=new double[11];

    for(int i=0;i<11;i++){
        cal[i]=callopition_montecarlo(sigma,r,T,s0[i],X);
        del[i]=delta(sigma,r,T,s0[i],X);
        the[i]=theta(sigma,r,T,s0[i],X);
        veg[i]=vega(sigma,r,T,s0[i],X);
        gam[i]=gamma(sigma,r,T,s0[i],X);
        rh[i]=rho(sigma,r,T,s0[i],X);
    }
    //output for graph
    ofstream file1;
    file1.open("C:\\Users\\Mia Dong\\Desktop\\Dong_Shuyu_Project3\\Shuyu_greeks.csv");
    for (int i = 0; i < 11; i++)
    {
        file1 <<s0[i]<<","<<cal[i]<<","<<del[i] << ","<<the[i]<<","<<veg[i]<<","<<gam[i]<<","<<rh[i]<< ",\n";
    }
    file1 << "\n";
    file1.close();

    //Q4
    //a)
    N=1000;
    double rho4=-0.6,r4=0.03,s04=48,v04=0.05,sigma4=0.42,alpha4=5.8,beta4=0.0625,T4=0.5,K4=50;
    double* HFT= new double[N];
    double*c1=new double[N];
    for(int i=0;i<N;i++){
        HFT[i]=Heston(sigma4,r4,T4,s04,v04,alpha4,beta4,rho4,0);
        c1[i]=exp(-r4*T4)*Callpayoff(HFT[i],K4);
    }
    double* HPT= new double[N];
    double*c2=new double[N];
    for(int i=0;i<N;i++){
        HPT[i]=Heston(sigma4,r4,T4,s04,v04,alpha4,beta4,rho4,1);
        c2[i]=exp(-r4*T4)*Callpayoff(HPT[i],K4);
    }
    double* HR= new double[N];
    double*c3=new double[N];
    for(int i=0;i<N;i++){
        HR[i]=Heston(sigma4,r4,T4,s04,v04,alpha4,beta4,rho4,2);
        c3[i]=exp(-r4*T4)*Callpayoff(HR[i],K4);
    }
    cout<<endl;
    cout<<"Q4: "<<endl;
    cout<<"Using the Full Truncation model, c1 is: "<<mean(c1,N)<<endl;
    cout<<"Using the Partial Truncation model, c2 is "<<mean(c2,N)<<endl;
    cout<<"Using the Reflection model, c3 is "<<mean(c3,N)<<endl;
    cout<<endl;

    //Q5
    //a)
    double*xa=get_uniform(100);
    double*ya=get_uniform(100);
    //b)
    double*xb=getHalton(100,2);
    double*yb=getHalton(100,7);
    //c)
    double*xc=getHalton(100,2);
    double*yc=getHalton(100,4);
    //d)output for graph
    ofstream file2;
    file2.open("C:\\Users\\Mia Dong\\Desktop\\Dong_Shuyu_Project3\\Shuyu_Haltons.csv");
    for (int i = 0; i < 100; i++)
    {
        file2 <<xa[i]<<","<<ya[i]<<","<<xb[i] << ","<<yb[i]<<","<<xc[i]<<","<<yc[i]<< ",\n";
    }
    file2 << "\n";
    file2.close();
    //e)
    N=10000;
    double*x1=getHalton(N,2);
    double*y1=getHalton(N,4);
    double*x2=getHalton(N,2);
    double*y2=getHalton(N,7);
    double*x3=getHalton(N,5);
    double*y3=getHalton(N,7);
    //(2,4)
    double*fun=new double[N];
    for(int i=0;i<N;i++){
        fun[i]=exp(-x1[i]*y1[i])*(sin(6*M_PI*x1[i])+cbrt(cos(2*M_PI*y1[i])));
    }
    cout<<"Q5: "<<endl;
    cout<<"Integral using bases(2,4) is: "<<mean(fun,N)<<endl;
    for(int i=0;i<N;i++){
        fun[i]=exp(-x2[i]*y2[i])*(sin(6*M_PI*x2[i])+cbrt(cos(2*M_PI*y2[i])));
    }
    cout<<"Integral using bases(2,7) is: "<<mean(fun,N)<<endl;
    for(int i=0;i<N;i++){
        fun[i]=exp(-x3[i]*y3[i])*(sin(6*M_PI*x3[i])+cbrt(cos(2*M_PI*y3[i])));
    }
    cout<<"Integral using bases(5,7) is: "<<mean(fun,N)<<endl;



    return 0;
}

    //c



























