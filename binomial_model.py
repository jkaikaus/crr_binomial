'''
Last Modified: 5/1/2017

Note - Python implementation of binomial model
'''

import numpy as np
import math
import sys
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from numpy import genfromtxt

#Python binomial implementation for runtime comparisons to C++
def Binomial(Option,K,T,S_not,sigma,r,q,N,Exercise):
    time_start = time.clock()

    #define all the variables
    delta = float(T)/float(N) #divide [0,T] into N equal intervals
    u = np.exp(sigma*np.sqrt(delta)) #calculate u and d based on delta
    d = np.exp(-1*sigma*np.sqrt(delta))
    p_star = float(np.exp(float(r-q)*delta)-d)/float(u-d)
    steps = np.linspace(0,N, num=N+1).astype(int)
    exponent = np.exp(-1*r*delta)
    q_star = 1-p_star

    #if option is European, so do not have to account for early exercise
    if Exercise == 'E':
        #backwards induction now
        if Option =='C':
            pay_off = (((u**steps)*(d**(N-steps))*S_not)-K).clip(min=0) #pricing of the last step
        elif Option == 'P':
            pay_off = (K-((u**steps)*(d**(N-steps))*S_not)).clip(min=0) #pricing of the last step
        else:
            return -1,-1 #error message, invalid option input
        for n in range(N-1,-1, -1):
            #temp = np.zeros(pay_off.shape)
            for j in range(0,n+1):
                pay_off[j] = exponent*(p_star*pay_off[j+1]+q_star*pay_off[j])

            #pay_off = np.copy(temp)
    #if option is American, now account for early exercise
    elif Exercise == 'A':
        if Option == 'C':
            pay_off = (((u**steps)*(d**(N-steps))*S_not)-K).clip(min=0)
            for n in range(N-1,-1,-1):
                for j in range(0,n+1):
                    S_j = (u**j)*(d**(n-j)*S_not)
                    pay_off[j] = max(max(0,S_j-K),exponent*(p_star*pay_off[j+1]+q_star*pay_off[j]))
        elif Option == 'P':
            pay_off = (K-((u**steps)*(d**(N-steps))*S_not)).clip(min=0)
            for n in range(N-1,-1,-1):
                for j in range(0,n+1):
                    S_j = (u**j)*(d**(n-j)*S_not)
                    pay_off[j] = max(max(0,K- S_j),exponent*(p_star*pay_off[j+1]+q_star*pay_off[j]))
        else:
            return -1, -1 #error message, invalid option input

    else:
        return -1,-1 #error message, option type

    return pay_off[0], time.clock() - time_start


#investigation of different scenarios
if __name__ == "__main__":
    #1 year euro call option with strike price of 100, current stock price of 100
    mydata = genfromtxt('problem2.csv', delimiter=",")
    black_scholes_price = black_scholes_check('C', 100, 1, 100,0.2,0.05,0.04)
    #plot results
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.title('Convergence of CRR Binomial Model')
    plt.plot(mydata[:,0],mydata[:,1], c='blue', label = 'CRR Binomial Model',linewidth=1.5)
    plt.axhline(y=black_scholes_price, linewidth=1.5, c='red', label = 'Black-Scholes Model')
    plt.xlabel('Time Steps')
    plt.ylabel('Option Price')
    plt.xlim((1,10000))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
          fancybox=True, shadow=True,ncol=2)
    #plt.show()
    fig.savefig('problem2.png')

    #american put option with strike price 100 dollars, finding the number of steps needed to get 10^-3 accuracy with time to maturity between 1 month to 1 year.
    put12monthdata = genfromtxt('problem3part1.csv', delimiter=",")
    put12monthdataq2 = genfromtxt('problem3part1q2.csv', delimiter=",")
    fig1 = plt.figure()
    ax = plt.subplot(111)
    plt.plot(put12monthdata[:,0], put12monthdata[:,1], c = 'b', label = 'q=0.00')
    plt.plot(put12monthdataq2[:,0], put12monthdataq2[:,1], c = 'r', label = 'q=0.04')
    plt.title('12 Month American Put: $S_0$ vs Option Price')
    plt.xlabel('Initial Stock Price')
    plt.ylabel('Option Price')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
          fancybox=True, shadow=True,ncol=2)
    #plt.show()
    fig1.savefig('problem3pt1.png')

    #early exercise boundary as function of time to maturity
    putearlyexercisedata = genfromtxt('problem3part2.csv', delimiter=",")
    putearlyexercisedataq2 = genfromtxt('problem3part2q2.csv', delimiter=",")
    fig2 = plt.figure()
    ax = plt.subplot(111)
    plt.plot(putearlyexercisedata[:,0], putearlyexercisedata[:,1], c='b',label = 'q=0.00')
    plt.plot(putearlyexercisedataq2[:,0], putearlyexercisedataq2[:,1], c='r', label = 'q=0.04')
    plt.title('American Put: Early Exercise Boundary')
    plt.xlabel('Time to Maturity (Months)')
    plt.ylabel('$S^*(i)$')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
          fancybox=True, shadow=True,ncol=2)
    #plt.show()
    fig2.savefig('problem3pt2.png')

    #american call option with strike price 100 dollars, finding the number of steps needed to get 10^-3 accuracy with time to maturity between 1 month to 1 year.
    call12monthdata = genfromtxt('problem4part1.csv', delimiter=",")
    call12monthdataq2 = genfromtxt('problem4part1q2.csv', delimiter=",")
    fig3 = plt.figure()
    ax = plt.subplot(111)
    plt.plot(call12monthdata[:,0], call12monthdata[:,1], c = 'b', label = 'q=0.04')
    plt.plot(call12monthdataq2[:,0], call12monthdataq2[:,1], c = 'r', label = 'q=0.08')
    plt.title('12 Month American Call: $S_0$ vs Option Price')
    plt.xlabel('Initial Stock Price')
    plt.ylabel('Option Price')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
          fancybox=True, shadow=True,ncol=2)
    #plt.show()
    fig3.savefig('problem4pt1.png')

    #critical stock price on the early exercise boundary
    callearlyexercisedata = genfromtxt('problem4part2.csv', delimiter=",")
    callearlyexercisedataq2 = genfromtxt('problem4part2q2.csv', delimiter=",")
    fig4 = plt.figure()
    ax = plt.subplot(111)
    plt.plot(callearlyexercisedata[:,0], callearlyexercisedata[:,1], c='b',label = 'q=0.04')
    plt.plot(callearlyexercisedataq2[:,0], callearlyexercisedataq2[:,1], c='r', label = 'q=0.08')
    plt.title('American Call: Early Exercise Boundary')
    plt.xlabel('Time to Maturity (Months)')
    plt.ylabel('$S^*(i)$')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
          fancybox=True, shadow=True,ncol=2)
    #plt.show()
    fig4.savefig('problem4pt2.png')
