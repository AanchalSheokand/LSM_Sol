# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 00:31:48 2025

@author: aanch
"""

def parameters():
    g=6.5
    pi=3.14
    c_1=4807.84
    h_l=1759728.23
    h_s=38072087.59
    m =-93795.1876
    lambda_1=13.49
    lambda_2=46.48
    fm=197.5
    return(g,pi,c_1,h_l,h_s,m,lambda_1,lambda_2,fm)

def mesonic_l(sigma_l,sigma_s):
    from math import sqrt
    from parameters__ import parameters
    (g,pi,c_1,h_l,h_s,m,lambda_1,lambda_2,fm)=parameters()
    u1=((2*lambda_1+lambda_2)/2)*(sigma_l**3)+(m-((sigma_s*c_1)/sqrt(2))+lambda_1*(sigma_s**2))*sigma_l-h_l
    return u1

def mesonics__(sigma_l,sigma_s):
    from parameters__ import parameters
    (g,pi,c_1,h_l,h_s,m,lambda_1,lambda_2,fm)=parameters()
    from math import sqrt
    
    u2=(lambda_1+lambda_2)*(sigma_s**3)+(m+lambda_1* (sigma_l**2))*sigma_s-h_s-((c_1*(sigma_l**2))/2*sqrt(2))
    return u2

def tempz_l(m_l,p_l,E_l):
    from math import log
    from parameters__ import parameters
    (g,pi,c,h_l,h_s,m,lambda_1,lambda_2,fm)=parameters()
    f1=((3/2)*((m_l*g)/pi**2)*((p_l*E_l)-(m_l**2)*log((p_l+E_l)/m_l)))
    return f1

def tempz_s(m_s,p_s,E_s):
    from math import log,sqrt
    from parameters__ import parameters
    (g,pi,c,h_l,h_s,m,lambda_1,lambda_2,fm)=parameters()
    
    f2=((3/(2*sqrt(2)))*((m_s*g)/pi**2)*((p_s*E_s)-(m_s**2)*log((p_s+E_s)/m_s)))
    return f2

from scipy.optimize import fsolve
from math import sqrt
import numpy as np 
from numpy import arctan,inf
from parameters__ import parameters
(g,pi,c_1,h_l,h_s,m,lambda_1,lambda_2,fm)=parameters()
sigmal = []
sigmas = []
rhobb = []
ml=[]
ms=[]
def fun(x,rhob):
    from scipy.integrate import quad
    from math import log,sqrt,pi
    from parameters__ import parameters
    (g,pi,c_1,h_l,h_s,m,lambda_1,lambda_2,fm)=parameters()

    
    
    sigma_l = x[0]
    sigma_s = x[1]
    
    

    s=0 # we give different values for strangeness fraction.
    
    
    rho_s = s*rhob
    rho_l = ((3-s)/2)*rhob
    
    
    p_l=(((pi**2)*rho_l))**(1/3)
    p_s = (((pi**2)*rho_s))**(1/3)      
    
    
    m_s = (g*sigma_s)/sqrt(2)
    m_l = (g*sigma_l)/2
    
    
    E_s = sqrt(p_s**2+m_s**2)
    E_l = sqrt(p_l**2+m_l**2)
    
    
    F1= mesonic_l(sigma_l, sigma_s)+tempz_l(m_l, p_l, E_l)
    F2= mesonics__(sigma_l, sigma_s)+tempz_s(m_s, p_s, E_s)
    
    
    F=[F1,F2]
    
    
    return F
for i in range(1,150 ):
    A = i*0.01
    rhob=A*(fm**3)
    s=0
    rho_s = s*A
    rho_l = ((3-s)/2)*A
    p_l=(((pi**2)*rho_l))**(1/3)
    p_s =(((pi**2)*rho_s))**(1/3) 

    
    x=fsolve(fun,(90,90),(rhob))

    rhobb.append((A/0.15)) # 0.15 here nuclear saturation density
    sigmal.append(x[0])
    sigmas.append(x[1])
    m_l=g*x[0]/2
    m_s=g*x[1]/sqrt(2)
    #print(m_l)
    #print(m_s)
    #print(A/0.15)
    
    ml.append(m_l)
    ms.append(m_s)
    print(sigmal)
    print(sigmas)
    print(rhobb)
    print(ml)
    print(ms)
   
    
   
    
import matplotlib   
import matplotlib.pyplot as plt
ax=plt.subplot()
matplotlib.rc('xtick',labelsize=18)
matplotlib.rc('ytick',labelsize=18)
plt.rc('legend',fontsize=14)


plt.xlabel(r"$ \rho /\rho_o$",fontsize=18)
#plt.ylabel("m(MeV)",fontsize=18)
plt.ylabel(r"$\sigma$(MeV)",fontsize=18)
ax.set_title(r'$\eta$=2',x=0.5,y=0.9,fontsize=18)
ax.set_xlim(0,10)
ax.set_ylim(0,110)
ax.set_ylim(0,500)
plt.plot(rhobb,ml,'--g')
plt.plot(rhobb,ms,'m')
#plt.plot(rhobb,sigmal,'--g')
#plt.plot(rhobb,sigmas,'m')
plt.legend([r"$\sigma_l{}$",r"$\sigma_s{}$"],loc="lower left")
#plt.legend(["$m_l$","$m_s$"],loc="lower left")

