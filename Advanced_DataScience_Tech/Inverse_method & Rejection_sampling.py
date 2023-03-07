from ast import arg
from cmath import exp
from lib2to3.pygram import Symbols
from pickletools import optimize
from random import uniform
from re import L
from statistics import NormalDist, mean
import numpy as np
import pandas as pd
import math
import scipy
import scipy.stats as st
import matplotlib.pyplot as plt
from scipy.special import gamma, factorial, psi, digamma
from scipy import optimize
from scipy.stats import kstest

#1
# Inverse function of "CDF_normal dist"
def normal_cdf_inverse(mu1, mu2, sigma1, sigma2):
    U=np.random.uniform(0,1,2)
    X = [0,0]
    X[0] = NormalDist(mu=mu1, sigma=sigma1).inv_cdf(U[0])
    X[1] = NormalDist(mu=mu2, sigma=sigma2).inv_cdf(U[1])
    return X
normal_cdf_inverse(0,0,1,1) # Generating two RVs from a Standard normal distribution
'''
We can input "mu1, mu2, sigma1, sigma2" which are the parameters of the Normal distribution that we want to generate sample from
'''

#2
#Drawing from Proposal function and Target function
#Using Inversion Method
def g(alpha):
    if alpha <= 0:
        print("alpha should be positive value")
    else:
        lambda_ = np.sqrt(2*alpha-1) if alpha>=1 else alpha
        mu = alpha**lambda_
        U=np.random.uniform(0,1,1)
        x = ((mu/(1-U))-mu)**(1/lambda_) # Inverse function of G, which is integral of g
        g= lambda_ * mu * x**(lambda_-1) / (mu+x**lambda_)**2
    return g,x # g(alpha)[0] is the RV generated from Proposal function

def f(alpha):
    if alpha <= 0:
        print("alpha should be positive value")
    else:
        x=g(alpha)[1][0]
        f= 1/gamma(alpha) * x**(alpha-1) * exp(-x)
    return f # f(alpha) is the RV generated from Target function

#Rejection Sampling
def rejection(alpha, n):
    if alpha <=0 or n<=0:
        print("alpha and 'n' should be positive value")
    y=[]
    while len(y) < n:
            x=g(alpha)[1][0]
            M=4*(alpha**alpha)*exp(-alpha) / (gamma(alpha)*np.sqrt(2*alpha-1)) if alpha>1 else 4*(alpha**alpha)*exp(-alpha) / (gamma(alpha)*alpha)
            #We already got M, at Hw1
            U=np.random.uniform(0,1,1)
            Ratio = f(alpha) / (M*g(alpha)[0][0])
            #g(alpha)[0][0]=RV generated from proposal function
            if U < Ratio:
                y.append(x)
                if len(y)==n:
                    break
    return y

def rejection_hist(alpha,n):
    plt.hist(rejection(alpha,n),range=(0,100),bins=50,alpha=0.5)
    plt.hist(np.random.gamma(alpha,1,n),range=(0,100),bins=50,alpha=0.5)
    plt.show()
rejection_hist(20,1000) #Example : alpha=20, n=1,000
'''
We can input "alpha" and "n=sample size" which are the parameters we want
We can also determine by 'ktest' whether samples can be said to have been generated from gamma distribution.
'''
kstest(rejection(20,1000), 'gamma', args=(20,0,1))
'''
KstestResult(statistic=0.08800554640731961, pvalue=3.452743409759758e-07)
pvalue is greater than 0.05, we can say samples are generated from gamma distribution
'''

#3
def gammaMLE(alpha, beta, n):
    if alpha <= 0 or beta <= 0 or n <= 0:
        print("Alpha, Beta and 'n' should be greater than zero")
    else:
        x=np.random.gamma(alpha,beta,n)
        x_bar = np.mean(x)
        
        def for_alpha(a):
            return -digamma(a)-np.log(beta)+np.mean(np.log(x))
        def for_beta(b):
            return -alpha/b + np.mean(x)/(b**2)
            
        alpha_hat = optimize.newton(for_alpha,1)
        beta_hat = optimize.newton(for_beta,0.5)
    
    return alpha_hat, beta_hat
gammaMLE(2,2,1000) #example : (1.9934150567308964, 2.0549092232225363)
'''
1)We can input "alpha, Beta, n" which are the parameters we want.
Then, samples(size n) are generated from the Gamma distribution.
2)We estimate alpha_hat and beta_hat assuming that we don't know from which distribution we produce samples.
We can see that alpht_hat and beta_hat are similar to the parameters we input beforehand.
'''

#4
def betaMLE(alpha, beta, n):
    if alpha <= 0 or beta <= 0 or n <= 0:
        print("Alpha, Beta and 'n' should be greater than zero")
    else:
        x=np.random.beta(alpha,beta,n)
        x_bar = np.mean(x)
        
        def for_alpha(a):
            return digamma(a+beta)-digamma(a)+np.mean(np.log(x))
        def for_beta(b):
            return digamma(alpha+b)-digamma(b)+np.mean(np.log(x))
            
        alpha_hat = optimize.newton(for_alpha,1)
        beta_hat = optimize.newton(for_beta,0.5)
    
    return alpha_hat, beta_hat
betaMLE(2,2,1000) #example : (2.0009085723245863, 2.0009085723244846)
'''
1)We can input "alpha, Beta, n" which are the parameters we want.
Then, samples(size n) are generated from the Beta distribution.
2)We estimate alpha_hat and beta_hat assuming that we don't know from which distribution we produce samples.
We can see that alpht_hat and beta_hat are similar to the parameters we input beforehand.
'''


#5
data = pd.read_csv("Lagrange.csv")
X=data["x"]
np.mean(X) # 0.029

def Lagrangian(X,iter):
    x_b = 0.5
    n = len(X)
    lambda_ = -1/(n*X)
    
    for i in range(iter):
        lambda_ = lambda_ - (x_b - np.sum(np.exp(-lambda_@X)*X))/(np.sum(np.exp(-lambda_@X)*X**2))
    w = np.exp(-lambda_@X)
    sample_mean = np.mean(X)
    weighted_mean = np.mean(w*100*X)
    print(f'simple mean is {sample_mean} and the weighted mean is {weighted_mean}')
    return w * 100

Lagrangian(data["x"].values, 100)
'''
simple mean is 0.029265651862 and the weighted mean is 0.05
1.7084874868248716
'''


#6
Raking = pd.read_csv("Raking.csv")
X1 = list(Raking["X1"])
X2 = list(Raking["X2"])

Raking_count = pd.DataFrame([(0,0,0),(0,0,0),(0,0,0)],index=("Male(1)","Female(2)","Sum"),columns=("Below(1)","Over(2)","Sum"))

for i in range(0,1000):
    if X1[i]==1 and X2[i]==1:
        Raking_count["Below(1)"][0] +=1
    elif X1[i]==1 and X2[i]==2:
        Raking_count["Over(2)"][0] +=1
    elif X1[i]==2 and X2[i]==1:
        Raking_count["Below(1)"][1] +=1
    elif X1[i]==2 and X2[i]==2:
        Raking_count["Over(2)"][1] +=1
Raking_count["Below(1)"][2]=Raking_count["Below(1)"][0]+Raking_count["Below(1)"][1]
Raking_count["Over(2)"][2]=Raking_count["Over(2)"][0]+Raking_count["Over(2)"][1]
Raking_count["Sum"][0]=Raking_count["Below(1)"][0]+Raking_count["Over(2)"][0]
Raking_count["Sum"][1]=Raking_count["Below(1)"][1]+Raking_count["Over(2)"][1]
Raking_count["Sum"][2]=Raking_count["Below(1)"][1]+Raking_count["Over(2)"][1]+Raking_count["Below(1)"][0]+Raking_count["Over(2)"][0]
Raking_count=Raking_count.astype(float)
Raking_count_original=Raking_count.copy()
'''
           Below(1)  Over(2)   Sum
Male(1)         376      216   592
Female(2)       135      273   408
Sum             511      489  1000
'''
#row-wise raking
def row_raking(DF_22):
    #row1
    DF_22["Below(1)"][0]=DF_22["Below(1)"][0]*500/DF_22["Sum"][0]
    DF_22["Over(2)"][0]=DF_22["Over(2)"][0]*500/DF_22["Sum"][0]
    DF_22["Sum"][0]=DF_22["Below(1)"][0]+DF_22["Over(2)"][0]
    #row2
    DF_22["Below(1)"][1]=DF_22["Below(1)"][1]*500/DF_22["Sum"][1]
    DF_22["Over(2)"][1]=DF_22["Over(2)"][1]*500/DF_22["Sum"][1]
    DF_22["Sum"][1]=DF_22["Below(1)"][1]+DF_22["Over(2)"][1]
    #col1 & col2
    DF_22["Below(1)"][2]=DF_22["Below(1)"][0]+DF_22["Below(1)"][1]
    DF_22["Over(2)"][2]=DF_22["Over(2)"][0]+DF_22["Over(2)"][1]
    return DF_22

def col_raking(DF_22):
    #col1
    DF_22["Below(1)"][0]=DF_22["Below(1)"][0]*300/DF_22["Below(1)"][2]
    DF_22["Below(1)"][1]=DF_22["Below(1)"][1]*300/DF_22["Below(1)"][2]
    DF_22["Below(1)"][2]=DF_22["Below(1)"][0]+DF_22["Below(1)"][1]
    #col2
    DF_22["Over(2)"][0]=DF_22["Over(2)"][0]*700/DF_22["Over(2)"][2]
    DF_22["Over(2)"][1]=DF_22["Over(2)"][1]*700/DF_22["Over(2)"][2]
    DF_22["Over(2)"][2]=DF_22["Over(2)"][0]+DF_22["Over(2)"][1]
    #row1 & row2
    DF_22["Sum"][0]=DF_22["Below(1)"][0]+DF_22["Over(2)"][0]
    DF_22["Sum"][1]=DF_22["Below(1)"][1]+DF_22["Over(2)"][1]
    return DF_22

row_raking_count=col_raking_count=0
while abs(Raking_count["Sum"][0]-500) > 10**(-6) or abs(Raking_count["Below(1)"][2] -300) > 10**(-6):
    row_raking(Raking_count)
    row_raking_count+=1
    col_raking(Raking_count)
    col_raking_count+=1

Raking_count
row_raking_count
col_raking_count
Raking_weights = Raking_count.div(Raking_count_original)
'''
>>> Raking_count
             Below(1)     Over(2)          Sum
Male(1)    212.961902  287.038097   499.999999
Female(2)   87.038098  412.961903   500.000001
Sum        300.000000  700.000000  1000.000000

#row_raking_count = 8
#col_raking_count = 8

>>> Raking_weights
           Below(1)   Over(2)       Sum
Male(1)    0.566388  1.328880  0.844595
Female(2)  0.644727  1.512681  1.225490
Sum        0.587084  1.431493  1.000000
->for example, 0.566388 = 212.961902 / 376
'''
