/*
Last Modified: 5/2/2017

Note - C++ (duh) implementation of binomial model
*/


#include <math.h>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <ctype.h>
#include <stdlib.h>
#include <string>
#include <time.h> //computation time calculations
#include <algorithm> //std::max
#include <fstream>  //outputing array data to txt file
#include <cmath>  //std::abs

struct RetVal {
  double price;
  double computational_time;
};

//check result with binomial implementation of Black-Scholes Model
RetVal Binomial(char Option, double K, double T, double S_not, double sigma, double r, double q, int N, char Exercise)
{
  RetVal price_time;
  clock_t t;
  t= clock();

  //variable definitions
  double delta = T/N;
  double u = exp(sigma*sqrt(delta));
  double d = exp(-1*sigma*sqrt(delta));
  double p_star = (exp((r-q)*delta)-d)/(u-d);
  double q_star = 1-p_star;
  double exponent = exp(-1*r*delta);
  std::vector<double> payoff;

  if (Exercise == 'E')
  {
    if (Option == 'C') {
      for (int i = 0; i<N+1; i++)
      {
        payoff.push_back(std::max(0.0,((pow(u,i)*pow(d,N-i))*S_not)-K));
      }
    } else if (Option == 'P') {
      for (int i = 0; i<N+1; i++)
      {
        payoff.push_back(std::max(0.0,K-((pow(u,i)*pow(d,N-i))*S_not)));
      }
    } else {
      std::cout << "Option must be either Call or Put!!!" << std::endl;
      price_time.price = -1;
      price_time.computational_time = (double)(clock()-t);
      return price_time;
    }

    for (int i = N-1; i >= 0; i--)
    {
      for (int j = 0; j<i+1;j++)
      {
        payoff.at(j) = exponent*(p_star*payoff.at(j+1)+q_star*payoff.at(j));
      }
    }

  } else if (Exercise == 'A'){
    if (Option == 'C') {
      for (int i = 0; i<N+1; i++)
      {
        payoff.push_back(std::max(0.0,((pow(u,i)*pow(d,N-i))*S_not)-K));
      }
      for (int i = N-1; i >= 0; i--)
      {
        for (int j = 0; j<i+1;j++)
        {
          double S_j = pow(u,j)*pow(d,i-j)*S_not;
          payoff.at(j) = std::max(std::max(0.0,S_j-K),exponent*(p_star*payoff.at(j+1)+q_star*payoff.at(j)));
        }
      }
    } else if (Option == 'P') {
      for (int i = 0; i<N+1; i++)
      {
        payoff.push_back(std::max(0.0,K-((pow(u,i)*pow(d,N-i))*S_not)));
      }
      for (int i = N-1; i >= 0; i--)
      {
        for (int j = 0; j<i+1;j++)
        {
          double S_j = pow(u,j)*pow(d,i-j)*S_not;
          payoff.at(j) = std::max(std::max(0.0,K-S_j),exponent*(p_star*payoff.at(j+1)+q_star*payoff.at(j)));
        }
      }
    } else {
      std::cout << "Option must be either Call or Put!!!" << std::endl;
      price_time.price = -1;
      price_time.computational_time = (double)(clock()-t);
      return price_time;
    }

  } else {
    std::cout << "Option must be either European or American!!!" << std::endl;
    price_time.price = -1;
    price_time.computational_time = (double)(clock()-t);
    return price_time;
  }

  price_time.price = payoff.front();
  price_time.computational_time = ((double)(clock()-t))/CLOCKS_PER_SEC;

  return price_time;
}

//investigation of different scenarios
int main(int argc, char *argv[])
{
  //1 year euro call option with strike price of 100, current stock price of 100
  double K = 100.0;
  double T = 1.0;
  double S = 100.0;
  int N = 10000;
  std::ofstream fout("../results1.csv");
  if (fout.is_open())
  {
    std::cout << "File open" << std::endl;
    //fout << 0 <<", "<< 0 << ", " << 0 << std::endl;
    for (int k= 100;k<=N;k+=100)
    {
      int idx = k/500;
      RetVal pt;
      pt = Binomial('C',K,T,S,0.2,0.05,0.04,k,'E');
      fout << k << ", "<< pt.price <<", " << pt.computational_time << std::endl;
      //std::cout<< "price for "<< k << ": "<<pt.price<< ", time: " << pt.computational_time <<  " seconds."<<std::endl;
    }
  } else {
    std::cout << "File could not be open" << std::endl;
  }

  //american put option with strike price 100 dollars, finding the number of steps needed to get 10^-3 accuracy with time to maturity between 1 month to 1 year.
  /*std::ofstream fout("../results2.csv");
  if (fout.is_open()){
    std::cout << "File open" << std::endl;
    double q1 = 0.0;
    double q2 = 0.04;
    double K = 100;
    double sigma = 0.2;
    double r = 0.05;
    double N = 2000; //number of steps needed to achieve 10^-3 accuracy
    double t12 = 1.0;
    for (double s = 75; s<=100; s++)
    {
      RetVal pt;
      pt = Binomial('P',K, t12, s, sigma, r, q2, N, 'A');
      fout << s << ", " << pt.price << ", " << pt.computational_time << std::endl;
    }
  } else {
    std::cout << "File could not be open" << std::endl;
  }

  //early exercise boundary as function of time to maturity
  double q1 = 0.0;
  double q2 = 0.04;
  double K = 100;
  double sigma = 0.2;
  double r = 0.05;
  double N = 2000;
  double tq1[13][4];
  for (int t = 0; t <= 12; t++)
  {
    double t_ = (double)t/12.0;
    tq1[t][0] = t;
    for (double s = 75; s<= 105; s+=0.01)
    {
      RetVal pt;
      pt = Binomial('P',K, t_, s, sigma, r, q2, N, 'A');
      double intrinsic = K - s;
      double abss = std::abs(intrinsic-pt.price);
      if (abss<0.005){
        tq1[t][1] = s;
        tq1[t][2] = pt.price;
        tq1[t][3] = pt.computational_time;
        //  std::cout << intrinsic << ", " << pt.price << std::endl;
        //  std::cout<< "Option price corresponding to stock price "<< s << " and time to maturity "<<t<< ":"<<pt.price<< ", time: " << pt.computational_time <<  " seconds. "<< abss<<std::endl;
      }
    }
    std::cout<< "Option price corresponding to stock price "<< tq1[t][1]<< " and time to maturity "<<tq1[t][0]<< ": "<<tq1[t][2]<< ", time: "<<tq1[t][3] << " seconds. "<<std::endl;

    std::cout<<"----------------------"<<std::endl;
    //break;
  }
  std::ofstream fout("../results3.csv");
  if (fout.is_open()){
    std::cout << "File open" << std::endl;
    for (int t = 0; t<=12; t++)
    {
      fout<<tq1[t][0]<<", "<< tq1[t][1]<<", " <<tq1[t][2]<<", "<<tq1[t][3]<<std::endl;
    }
  }

  //american call option with strike price 100 dollars, finding the number of steps needed to get 10^-3 accuracy with time to maturity between 1 month to 1 year.
   std::ofstream fout("../results4.csv");
   if (fout.is_open()){
    std::cout << "File open" << std::endl;
    double q1 = 0.04;
    double q2 = 0.08;
    double K = 100;
    double sigma = 0.2;
    double r = 0.05;
    double N = 2000; //number of steps needed to achieve 10^-3 accuracy
    double t12 = 1.0;
    for (double s = 100; s<=140; s++)
    {
      RetVal pt;
      pt = Binomial('C',K, t12, s, sigma, r, q1, N, 'A');
      fout << s << ", " << pt.price << ", " << pt.computational_time << std::endl;
    }
  } else {
    std::cout << "File could not be open" << std::endl;
  }

  //critical stock price on the early exercise boundary
  double q1 = 0.04;
  double q2 = 0.08;
  double K = 100;
  double sigma = 0.2;
  double r = 0.05;
  double N = 2000;
  double tq1[13][4];
  for (int t = 0; t <= 12; t++)
  {
    double t_ = (double)t/12.0;
    tq1[t][0] = t;
    double start = 95;
    for (double s = start; s<= 170; s+=0.01)
    {
      RetVal pt;
      pt = Binomial('C',K, t_, s, sigma, r, q2, N, 'A');
      double intrinsic = s-K;
      double abss = std::abs(intrinsic-pt.price);
      if (abss<0.005){
        tq1[t][1] = s;
        tq1[t][2] = pt.price;
        tq1[t][3] = pt.computational_time;
        start = s;
          //std::cout << intrinsic << ", " << pt.price << std::endl;
          std::cout<< "Option price corresponding to stock price "<< s << " and time to maturity "<<t<< ":"<<pt.price<< ", time: " << pt.computational_time <<  " seconds. "<< abss<<std::endl;
        break;
      }
    }
    //break;
    //std::cout<< "Option price corresponding to stock price "<< tq1[t-1][1]<< " and time to maturity "<<tq1[t-1][0]<< ":"<<tq1[t-1][2]<< ", time: " << " seconds. "<<std::endl;

    std::cout<<"----------------------"<<std::endl;
    //break;
  }
  std::ofstream fout("../results5.csv");
  if (fout.is_open()){
    std::cout << "File open" << std::endl;
    for (int t = 0; t<=12; t++)
    {
      fout<<tq1[t][0]<<", "<< tq1[t][1]<<", " <<tq1[t][2]<<", "<<tq1[t][3]<<std::endl;
    }
  }

  //miscellaneous
  //10^-3 accuracy for steps for american calls and puts
  // double K = 100.0;
  // double T = 1.0;
  // double S = 100.0;
  // int N = 10000;
  // std::ofstream fout("../results6.csv");
  // if (fout.is_open())
  // {
  //   std::cout << "File open" << std::endl;
  //   //fout << 0 <<", "<< 0 << ", " << 0 << std::endl;
  //   for (int k= 500;k<=N;k+=500)
  //   {
  //     int idx = k/500;
  //     RetVal pt;
  //     pt = Binomial('P',K,T,S,0.2,0.05,0.04,k,'E');
  //     fout << k << ", "<< pt.price <<", " << pt.computational_time << std::endl;
  //     //std::cout<< "price for "<< k << ": "<<pt.price<< ", time: " << pt.computational_time <<  " seconds."<<std::endl;
  //   }
  // } else {
  //   std::cout << "File could not be open" << std::endl;
  // }
  */
}
