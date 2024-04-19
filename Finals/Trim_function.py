#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 14:08:15 2024

@author: nathanrawiri
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from sympy import symbols, solve
from numpy import cos, sin


def D_(C_D, rho, V, S):
    return C_D * 0.5 * rho * V ** 2 * S


def C_T_(T, rho, omega, R):
    return T / (rho * (omega * R) ** 2 * math.pi * R ** 2)


def T_(W, D):
    return math.sqrt(W ** 2 + D ** 2)


"""initilise variables"""

rho = 1.225  # [kg / m^3]
g = 9.81  # [m/s^2]
m_emp = 1051  # [kg]
m_fuel = 405  # [kg]
m = m_emp + m_fuel  # [kg]
W = m * g  # [N]
CDfus = 0.03
Sfus = 38.585
R = 5.345  # [m]
omega = 213 / R  # rad/s
Cla = 5.7
n_blades = 3
blade_chord = 0.3  # [m]
sigma = n_blades * blade_chord / (math.pi * R)


def getTrimAngles(Vh):
    vi = symbols('vi')
    """initialise variables"""
    meu = Vh / (omega * R)
    D = D_(CDfus, rho, Vh, Sfus)
    T = T_(W, D)
    C_T = C_T_(T, rho, omega, R)

    """The below 3 lines is what causes the program to run so slowly. sp solve"""
    error = 10 ** (-18)
    lambdai0 = 0.1
    lambdai = 0
    iter = 1000
    Ct = C_T
    for i in range(iter):
        V1 = (Vh / (omega * R) * cos(D / W))
        V2 = (Vh / (omega * R) * sin(D / W) + lambdai0)
        lambdai = (Ct / 2) / np.sqrt(V1 ** 2 + V2 ** 2)
        if (np.abs((lambdai - lambdai0) / lambdai)) < error:
            break
        else:
            lambdai0 = lambdai
    vibar = lambdai
    """Above is the slow part"""

    """Do matrix manipulation"""
    """Take matrix equation as A(2x2).B(1x2) = C(1x2). Rearrange for B = A-1 . C"""
    A = [[1 + (3 / 2) * meu ** 2, -8 / 3 * meu],
         [-meu, 2 / 3 + meu ** 2]]
    C = [[-2 * meu ** 2 * D / W - 2 * meu * vibar],
         [4 / sigma * C_T / Cla + meu * D / W + vibar]]
    A_inv = np.linalg.inv(A)
    B = np.dot(A_inv, C)

    """Get collective and cyclic angles"""
    cyclic = B[0][0]
    collect = B[1][0]
    # print(f"cyclic = {cyclic *180/np.pi}")
    # print(f"collect = {collect *180/np.pi}")
    # print(f"Pitch = {theta_f}")
    return cyclic, collect, lambdai
