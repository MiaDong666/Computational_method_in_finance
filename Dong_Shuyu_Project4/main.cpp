#include <iostream>
#include <vector>
#include "head.h"
#include <sstream>
#include <fstream>
#include<cmath>

using namespace std;

int main() {
    //Q1
    double r=0.05,t=0.5,sigma=0.24,s0=32,K=30;
    double* C1=new double[7];
    double* C2=new double[7];
    double* C3=new double[7];
    double* C4=new double[7];
    vector<int> steps={10,20,40,80,100,200,500};
    for(int i=0;i<7;i++){
        C1[i]=Binomial_Europeanoption(steps[i],sigma,r,t,s0,K,0);
        C2[i]=Binomial_Europeanoption(steps[i],sigma,r,t,s0,K,1);
        C3[i]=Binomial_Europeanoption(steps[i],sigma,r,t,s0,K,2);
        C4[i]=Binomial_Europeanoption(steps[i],sigma,r,t,s0,K,3);
    }

    ofstream file1;
    file1.open("C:\\Users\\Mia Dong\\Desktop\\HW4\\Shuyu_Q1.csv");
    for (int i = 0; i < 7; i++)
    {
        file1 <<C1[i]<<","<<C2[i]<<","<<C3[i] << ","<<C4[i]<< ",\n";
    }
    file1 << "\n";
    file1.close();


    //Q2
    //get the historical sd
    ifstream  data("C:\\Users\\Mia Dong\\Desktop\\HW4\\goog.csv");
    string line;
    int counter = 0;
    double *ret = new double[59];

    while(getline(data,line))
    {
        long double temp = stold(line);
        ret[counter] = temp;
        counter ++;
    }
    sigma=sqrt(12*var(ret,59));
    //GOOG PRICE  AS OF 02/05/2019
    r=0.02,s0=1146,K=1260,t=1;
    cout<<"Q2: "<<endl;
    cout<<Binomial_Europeanoption(500,sigma,r,t,s0,K,3)<<endl;
    //the call option expires 2020 at K=1260 is $56.9
    //to calculate implied volatility, we simulate by step.
    for(int i=0;i<1000;i++){
        if(Binomial_Europeanoption(500,sigma,r,t,s0,K,3)<56.9){
            sigma=sigma+0.001;
        }else{
            break;
        }
    }
    cout<<"implied volatility is: "<<sigma<<endl;
    cout<<endl;

    //Q3
    s0=49,K=50,r=0.03,sigma=0.2,t=0.3846;
    double *del1= new double[31];
    double *the=new double[31];
    double *gam=new double[31];
    double *veg=new double[31];
    double*rh=new double[31];

    double S=20;
    double epsilon=0.01;
    for(int i=0;i<31;i++) {
        del1[i] = (Binomial_Europeanoption(500,sigma, r, t, S+epsilon, K,3)-
                Binomial_Europeanoption(500,sigma, r, t, S-epsilon, K,3))/(2*epsilon);
        the[i] = -(Binomial_Europeanoption(500,sigma, r, t+epsilon, S, K,3)-
                  Binomial_Europeanoption(500,sigma, r, t-epsilon, S, K,3))/(2*epsilon);
        gam[i] = (Binomial_Europeanoption(500,sigma, r, t, S+1, K,3)+
                  Binomial_Europeanoption(500,sigma, r, t, S-1, K,3)-
                  2*Binomial_Europeanoption(500,sigma, r, t, S, K,3))/(1*1);
        veg[i] = (Binomial_Europeanoption(500,sigma+epsilon, r, t, S, K,3)-
                  Binomial_Europeanoption(500,sigma-epsilon, r, t, S, K,3))/(2*epsilon*100);
        rh[i] = (Binomial_Europeanoption(500,sigma, r+epsilon, t, S, K,3)-
                 Binomial_Europeanoption(500,sigma, r-epsilon, t, S, K,3))/(2*epsilon*100);
        S = S + 2;
    }
    ofstream file2;
    file2.open("C:\\Users\\Mia Dong\\Desktop\\HW4\\Shuyu_Q3_1.csv");
    for (int i = 0; i < 31; i++)
    {
        file2 <<del1[i]<<","<<the[i]<<","<<gam[i] << ","<<veg[i]<<","<<rh[i]<< ",\n";
    }
    file2 << "\n";
    file2.close();

    double * del2=new double[39];
    double t0=0;
    for(int i=0;i<39;i++) {
        del2[i]=(Binomial_Europeanoption(500,sigma, r, t0, s0+1, K,3)-
                 Binomial_Europeanoption(500,sigma, r, t0, s0-1, K,3))/(2*1);
        t0 = t0+0.01;
    }
    ofstream file3;
    file3.open("C:\\Users\\Mia Dong\\Desktop\\HW4\\Shuyu_Q3_2.csv");
    for (int i = 0; i < 39; i++)
    {
        file3 <<del2[i]<< ",\n";
    }
    file3 << "\n";
    file3.close();

    //Q4
    r=0.05,sigma=0.3,K=100,t=1;
    S=80;
    double *put1= new double[11];
    double *put2=new double[11];
    for(int i=0;i<11;i++){
        put1[i]=Binomial_put(500,sigma, r, t, S, K,0);
        put2[i]=Binomial_put(500,sigma, r, t, S, K,1);
        S=S+4;
    }
    ofstream file4;
    file4.open("C:\\Users\\Mia Dong\\Desktop\\HW4\\Shuyu_Q4.csv");
    for (int i = 0; i < 11; i++)
    {
        file4 <<put1[i]<<","<<put2[i]<< ",\n";
    }
    file4 << "\n";
    file4.close();;
















    return 0;
}