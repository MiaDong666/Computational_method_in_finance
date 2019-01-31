//computational methods
//HW1

#include "head.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<random>
#include<stdio.h>
#include <algorithm>
#include <chrono>
using namespace std;
using namespace std::chrono;


int main()
{
	//Q1 a) LGM
	cout<<"Q1 answer:"<<endl;
	int n=10000;

	double * result_LGM = get_uniform(n);
	cout<<"mean_LGM is: "<< mean(result_LGM, n) <<endl;
	cout << "sd_LGM is: " << sqrt(var(result_LGM, n)) << endl;

	//Q1 b) built-in function
	double* result_UN = get_uniform_UN(n);
	cout << "mean_UN is: " << mean(result_UN, n) << endl;
	cout << "sd_UN is: " << sqrt(var(result_UN, n)) << endl;

	//output for graph
	ofstream file1;
	file1.open("C:\\Users\\Mia Dong\\Desktop\\computational_hw1\\Shuyu_Uniform.csv");
	for (int i = 0; i < n; i++)
	{
		file1 << result_LGM[i] << "," << result_UN[i] <<  ",\n";
	}
	file1 << "\n";
	file1.close();

	//Q2 a)
	cout << "Q2 answer:" << endl;
	double * result_X=get_uniform(n);
	for (int i=0;i<n;i++) 
	{
		if (result_X[i]>=0 && result_X[i]<0.3) {
			result_X[i] = -1;
		}
		else if(result_X[i]>=0.3 && result_X[i]<0.65)
		{
			result_X[i] = 0;
		}
		else if (result_X[i]>=0.65 && result_X[i]<0.85) 
		{
			result_X[i] = 1;		
		}
		else 
		{
			result_X[i] = 2;
		}
	}
	//Q2 b)
	cout << "mean_X is: "<< mean(result_X, n) << endl;
	cout << "sd_X is: " << sqrt(var(result_X, n)) << endl;

	//output for graph
	ofstream file2;
	file2.open("C:\\Users\\Mia Dong\\Desktop\\computational_hw1\\Shuyu_X.csv");
	for (int i = 0; i < n; i++)
	{
		file2 << result_X[i] << ",\n";
	}
	file2 << "\n";
	file2.close();


	//Q3 a)
	cout << "Q3 answer:" << endl;
	double* result_bi;
	double p = 0.64;
	int k = 44;
	n = 1000;
	result_bi = Binom(p,k,n);//get binomial RV
	
	//Q3 b)calculate probability
	cout<<"probability X>=40 is: "<<probbig(result_bi,n,40)<<endl;
	
	//output for graph
	ofstream file3;
	file3.open("C:\\Users\\Mia Dong\\Desktop\\computational_hw1\\Shuyu_binomial.csv");
	for (int i = 0; i < n; i++)
	{
		file3 << result_bi[i] << ",\n";
	}
	file3 << "\n";
	file3.close();

	//Q4 a)
	cout << "Q4 answer:" << endl;
	double *result_exp;
	double lambda = 1.5;
	n = 10000;
	result_exp = Exp(lambda, n);

	//Q4 b)
	cout << "Probability that X>=1 is: "<<probbig(result_exp,n,1) <<endl;
	cout << "Probability that X>=4 is: "<<probbig(result_exp,n,4) <<endl;

	//Q4 c)
	cout << "mean_exp is: " << mean(result_exp,n)<<endl;
	cout << "sd_exp is: " << sqrt(var(result_exp, n)) << endl;
	   
	//output for graph
	ofstream file4;
	file4.open("C:\\Users\\Mia Dong\\Desktop\\computational_hw1\\Shuyu_exp.csv");
	for(int i = 0; i < n; i++)
	{
		file4 <<result_exp[i] << ",\n";
	}


	file4.close();
	
	//Q5 a)
	cout<<"Q5 answer:"<<endl;
	n = 5000;
	double * result_U;
	result_U = get_uniform(n);

	//Q5 b)
	double*result_N_BM;
	result_N_BM=normal_BM(n);

	//Q5 c)
	cout<<"the mean of Box-Muller Method normal is: "<<mean(result_N_BM,n)<<endl;
	cout<<"the standard deviation of Box-Muller normal is: "<<sqrt(var(result_N_BM,n))<<endl;

    //Q5 d)
    double *result_N_PM;
    result_N_PM=normal_PM(n);

    //Q5 e)
    cout<<"the mean of Polar Marsaglia normal is: "<<mean(result_N_PM,n)<<endl;
    cout<<"the standard deviation of Polar Marsaglia normal is: "<<sqrt(var(result_N_PM,n))<<endl;

    //Q5 f)
    double*test;
    int N=1000;
    int c_BM=0;
    for (int i=0;i<N;i++)
    {
    	//time for both normal generating functions
    	//BM
    	double t_BM,t_PM;
        auto start1=high_resolution_clock::now();
        test=normal_BM(n);
        auto stop1=high_resolution_clock::now();
        auto duration1=duration_cast<microseconds>(stop1-start1);
        t_BM=duration1.count();
        //PM
		auto start2=high_resolution_clock::now();
		test=normal_PM(n);
		auto stop2=high_resolution_clock::now();
		auto duration2=duration_cast<microseconds>(stop2-start2);
		t_PM=duration2.count();

		//compare which one runs faster
		if(t_BM<t_PM)
		{
			c_BM=c_BM+1;
		}
    }
    cout<<"in "<<N<<" times of test, Box-Muller runs faster than Polar Marsaglia for "<<c_BM<<" times"<<endl;



    //output for graph
    ofstream file5;
    file5.open("C:\\Users\\Mia Dong\\Desktop\\computational_hw1\\Shuyu_Normal.csv");
    for (int i = 0; i < n; i++)
    {
        file5 << result_N_BM[i] << "," << result_N_PM[i] <<  ",\n";
    }
    file5 << "\n";
    file5.close();


	
	return 0;

}



