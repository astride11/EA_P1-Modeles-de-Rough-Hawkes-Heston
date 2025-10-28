#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 2020

@author: Sergio Pulido
"""
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import scipy.integrate as integrate
from scipy.special import gamma

# Characteristic function of the Lifted Heston model see Slides 85-87
def Ch_Lifted_Heston(omega,S0,T,rho,lamb,theta,nu,V0,N,rN,alpha,M):
    # omega = argument of the ch. function
    # S0 = Initial price
    # rho,lamb,theta,nu,V0 = parameters Lifted Heston
    # N = number of factors in the model
    # rN = constant used to define weights and mean-reversions
    # alpha = H+1/2 where H is the Hurst index
    # T = maturity
    # M = number of steps in the time discretization to calculate ch. function

    # to make sure we calculate ch. function and not moment gen. function
    i=complex(0,1)
    omega=i*omega
    
    # Definition of weights and mean reversions in the approximation
    h=np.linspace(0,N-1,N)
    rpowerN=np.power(rN,h-N/2) 
    # weights
    c=(rN**(1-alpha)-1)*(rpowerN**(1-alpha))/(gamma(alpha)*gamma(2-alpha))
    # mean reversions 
    gammas=((1-alpha)/(2-alpha))*((rN**(2-alpha)-1)/(rN**(1-alpha)-1))*rpowerN
    
    # Definition of the initial curve
    g = lambda t: V0+lamb*theta*np.dot(c/gammas,1-np.exp(-t*gammas))
    
    
    # Time steps for the approximation of psi         
    delta = T/M
    t=np.linspace(0,M,M+1)
    t = t * delta
    
    # Function F
    F = lambda u,v : 0.5*(u**2-u)+(rho*nu*u-lamb)*v+.5*nu**2*v**2
    
    
    # Iteration for approximation of psi - see Slide 87
    psi=np.zeros((M+1,N),dtype=complex)
    
    for k in range (1,M+1):
        psi[k,:] = (np.ones(N)/(1+delta*gammas))*(psi[k-1,:]+delta*F(omega,np.dot(c,psi[k-1,:]))*np.ones(N))
        
    
    # Invert g_0 to calculate phi - see Slide 87
    g_0=np.zeros((1,M+1))
    
    for k in range(1,M+2):
        g_0[0,k-1]=g(T-t[k-1])
    
    
    Y=np.zeros((1,M+1),dtype=complex)
    phi=0
    
    Y=F(omega,np.dot(c,psi.transpose()))*g_0
   
    
    # Trapezoid rule to calculate phi
    weights=np.ones(M+1)*delta
    weights[0]=delta/2
    weights[M]=delta/2
    phi=np.dot(weights,Y.transpose())
    
    phi=np.exp(omega*np.log(S0)+phi)
    
    return phi
    
    
# Computation of the price of Call option in the Lifted Heston model using Lewis formula
def Call_Price_Lifted_Heston(S0,K,T,rint,rho,lamb,theta,nu,V0,N,rN,alpha,M,alpha2,L):
    # S0 = Initial price
    # K,T = strike and maturity of the call option
    # rint = risk free rate
    # rho,lamb,theta,nu,V0 = parameters Lifted Heston
    # N = number of factors in the model
    # rN = constant used to define weights and mean-reversions
    # alpha = H+1/2 where H is the Hurst index
    # M = number of steps in the time discretization to calculate ch. function
    # alpha2 = damping factor in the Carr-Madan formula
    # L = truncation bound for the integral

    # Lifted Heston characteristic function 
    i=complex(0,1)
    phi=lambda omega:Ch_Lifted_Heston(omega,S0,T,rho,lamb,theta,nu,V0,N,rN,alpha,M)
    
    
    #Integrand in Carr-Madan formula - see notes on the Fourier transform fmla. (4.36)-(4.37)
    integrand=lambda omega:np.real((phi(omega-i*(alpha2+1))/(alpha2**2+alpha2-omega**2+i*(2*alpha2+1)*omega))*np.exp(-i*np.log(K)*omega))                 
    
    # Pricing formula - Carr-Madan - see notes on the Fourier transform fmla. (4.36)-(4.37)
    I = integrate.quad(integrand,0,L)
    P=(np.exp(-rint*T-alpha2*np.log(K))/np.pi)*np.asarray(I)[0] 
    
    return P


from scipy import optimize,stats

#Compute imlied volatility
def BS_pricer(K,S0, T_annual, r, sigma) :
    d1 = 1/(sigma*np.sqrt(T_annual))*(np.log(S0/K)+(r+1/2*sigma**2)*T_annual)
    d2 = d1 - sigma*np.sqrt(T_annual)
    Call_price = S0*stats.norm.cdf(d1) - K*np.exp(-r*T_annual)*stats.norm.cdf(d2)
    return Call_price


def imp_vol(K, S0, r, T_annual, Call_price) :
    def sq_cost(sigma) :
        cost = (Call_price- BS_pricer(K, S0,  T_annual, r, sigma[0]))**2
        return cost
    res = optimize.differential_evolution(sq_cost, bounds =[(0,1)], maxiter=100)
    return float(res['x'])
    
S0=1;rint=0;rho=-.7;lamb=.3;theta=.02;nu=.3;V0=.02;N=20
rN=1+10*(N**(-.9));alpha=.6;M=100;L=1000;alpha2=1

T=1/26

log_K=np.linspace(-.15,.05,20)
Nk=np.size(log_K)

IV=np.zeros(Nk)
P=np.zeros(Nk)

for i in range(1,Nk+1):
    P[i-1]=Call_Price_Lifted_Heston(1,np.exp(log_K[i-1]),T,rint,rho,lamb,theta,nu,V0,N,rN,alpha,M,alpha2,L)
    IV[i-1]=imp_vol(np.exp(log_K[i-1]), S0, rint, T, P[i-1])
    
import matplotlib.pyplot as plt

plt.plot(log_K,IV)

