'''
Last Modified: 4/20/2017

Note - Black-Scholes model implementation. Used to check for convergence on c++ and python implementations of binomial model
'''


import numpy as np
import math
import sys
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from numpy import genfromtxt

def black_scholes_check(Option,K,T,S_not,sigma,r,q):
    d1 = float(np.log(float(S_not)/float(K)) + (r-q+.5*sigma**2)*T)/float(sigma*np.sqrt(T))
    print 'd1: ', d1
    d2 = d1 - sigma*np.sqrt(T)
    print 'd2: ', d2
    price = 0
    if Option == 'C':
        price = S_not*np.exp(-1*q*T)*norm.cdf(d1) - K*np.exp(-r*T)*norm.cdf(d2)
    else:
        print norm.cdf(-1*d1)
        print norm.cdf(-1*d2)
        price = -1*S_not*np.exp(-1*q*T)*norm.cdf(-1*d1) + K*np.exp(-r*T)*norm.cdf(-1*d2)
    return price


if __name__ == "__main__":
    print black_scholes_check('P',1.58,1,1.56,.2,.02,.03)
