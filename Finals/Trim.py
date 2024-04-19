#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 14:08:15 2024

Description: 
This program creates a plots 
collective and cyclic pitch angles in trim 
at various velocities from hover to max V
for our chosen helicopter.

@authors: Nathan Rawiri and Matteus Ramos
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from sympy import symbols, solve
import sympy as sp


def D_(C_D, rho, V, S):
    return C_D*0.5*rho*V**2*S

def C_T_(T, rho, omega, R):
    return T / (rho*(omega*R)**2*math.pi*R**2)
    
def T_(W, D):
    return math.sqrt(W**2 + D**2)


"""initilise variables"""
rho = 1.225                                         # [kg / m^3]
g = 9.81                                            # [m/s^2]
m_emp = 1051                                        # [kg]
m_fuel = 405                                        # [kg]
m = m_emp + m_fuel                                  # [kg]
W = m*g                                             # [N]
CDfus = 0.03 
Sfus = 38.585                                       # [m^2]  
R = 5.345                                           # [m]
omega = 213/R                                       # rad/s
Cla = 5.7
n_blades = 3
blade_chord = 0.3                                   # [m]
sigma = n_blades * blade_chord / (math.pi * R)



"""For a given velocity, this function returns the cyclic and collective pitch for trim"""
def getTrimAngles(Vh):
    
    # Initialise variables
    meu = Vh / (omega * R)
    D = D_(CDfus, rho, Vh, Sfus)
    T = T_(W, D)
    C_T = C_T_(T, rho, omega, R)
    
    # Get lambda_i
    lambi = symbols('lambi')
    """The below 3 lines is what causes the program to run so slowly. sp solve"""
    lambi_equation = (2*lambi*sp.sqrt((Vh/(omega*R)*sp.cos(D/W))**2+(Vh/(omega*R)*sp.sin(D/W)+lambi)**2)) - C_T
    lambda_i = solve(lambi_equation, lambi) #"""This line is what takes a long time"""
    vibar = lambda_i[0]
    
    
    """Do matrix manipulation"""
    # - Take matrix equation as A(2x2).B(1x2) = C(1x2). Rearrange for B = A-1 . C
    A = [[1+(3/2)*meu**2,-8/3*meu],
        [-meu, 2/3 + meu**2]]
    C = [[-2*meu**2*D/W-2*meu*vibar],
          [4/sigma*C_T/Cla+meu*D/W+vibar]]
    A_inv = np.linalg.inv(A)
    B = np.dot(A_inv, C)
    
    """Get collective and cyclic angles"""
    cyclic = B[0][0]
    collect = B[1][0]
    #print(f"cyclic = {cyclic *180/np.pi}")
    #print(f"collect = {collect *180/np.pi}")
    return cyclic, collect

def main():
    cyclics = []
    collects = []
    
    num_iterations = 40
    Vhs = np.linspace(0.1, 80, num_iterations)
    for Vh in Vhs:
        print(f"V = {Vh:.2f}")
        cyclic, collect = getTrimAngles(Vh)
        cyclics.append(cyclic * 180/np.pi)
        collects.append(collect *180/np.pi)
    
    plt.figure()
    plt.xlabel('V (m/s)')
    plt.ylabel('θ')
    plt.plot(Vhs, collects, label="collective (θ0)", color="b")
    plt.plot(Vhs, cyclics, label="cyclic (θc)", color="c")
    plt.legend()
    plt.grid()
    plt.show()

main()
print(f"Trim for 90 knots:  Cyclic = {getTrimAngles(46.3)[0]*180/np.pi:.2f} deg,\
      Collective =  {getTrimAngles(46.3)[1]*180/np.pi:.2f} deg")
